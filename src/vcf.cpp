#include <string>
#include <unordered_map>
#include <vector>

#include "htslib/vcf.h"

#include "vcf.h"
#include "print.h"

/******************************************************************************/

variantCalls::variantCalls() { this->offs.resize(2); };

void variantCalls::add_gap(int g) { this->gaps.push_back(g); }

void variantCalls::add_var(int pos, int rlen, int off0, int off1, uint8_t hap, 
        uint8_t type, std::string ref, std::string alt, float gq, float vq) {
    this->poss.push_back(pos);
    this->rlens.push_back(rlen);
    this->offs[0].push_back(off0);
    this->offs[1].push_back(off1);
    this->haps.push_back(hap);
    this->types.push_back(type);
    this->refs.push_back(ref);
    this->alts.push_back(alt);
    this->gt_quals.push_back(gq);
    this->var_quals.push_back(vq);
}

/******************************************************************************/

vcfData::vcfData(htsFile* vcf) : hapcalls(2) {

    // counters
    int nseq   = 0;                     // number of sequences
    std::vector< std::vector<int> > 
        ntypes(2, std::vector<int>(type_strs.size(), 0));
    int n      = 0;                     // total number of records in file
    int npass  = 0;                     // records PASSing all filters

    // data
    std::vector<int> prev_end = {-g.gap*2, -g.gap*2};
    std::unordered_map<int, bool> prev_rids;
    int prev_rid = -1;
    std::string seq;

    // need two versions, since grouping is different if hap only vs both
    int var_idx = 0;                       // variant indices
    std::vector<int> hap_var_idx = {0, 0}; // indices per hap/contig
    std::vector<int> offs = {0, 0};        // offsets in group after INDELs
    std::vector<int> hap_offs = {0, 0};    // offsets in group after INDELs

    // quality data for each call
    int ngq_arr = 0;
    int ngq     = 0;
    int * gq    = (int*) malloc(sizeof(int));
    float * fgq = (float*) malloc(sizeof(float));
    bool int_qual = true;

    // genotype data for each call
    int ngt_arr   = 0;
    bool gq_warn  = false;
    int ngt       = 0;
    std::vector<int> ngts(gt_strs.size(), 0);
    int * gt      = NULL;
    
    // read header
    bcf1_t * rec  = NULL;
    bcf_hdr_t *hdr = bcf_hdr_read(vcf);
    int pass_filter_id = 0;
    bool pass = false;
    bool pass_found = false;
    for(int i = 0; i < hdr->nhrec; i++) {

        // search all FILTER lines
        if (hdr->hrec[i]->type == BCF_HL_FLT) {

            // select PASS filter
            bool is_pass_filter = false;
            for(int j = 0; j < hdr->hrec[i]->nkeys; j++)
                if (hdr->hrec[i]->keys[j] == std::string("ID") && 
                        hdr->hrec[i]->vals[j] == std::string("PASS"))
                    is_pass_filter = true;

            // save PASS filter index to keep only passing reads
            if (is_pass_filter)
                for(int j = 0; j < hdr->hrec[i]->nkeys; j++)
                    if (hdr->hrec[i]->keys[j] == std::string("IDX")) {
                        pass_filter_id = std::stoi(hdr->hrec[i]->vals[j]);
                        pass_found = true;
                    }
        }
    }
    if (!pass_found) {
        ERROR("failed to find PASS FILTER in VCF");
    }

    // report names of all the sequences in the VCF file
    const char **seqnames = NULL;
    seqnames = bcf_hdr_seqnames(hdr, &nseq);
    if (seqnames == NULL) {
        ERROR("failed to read VCF header");
        goto error1;
    }
    for(int i = 0; i < nseq; i++) {
        this->hapcalls[0][seqnames[i]] = variantCalls();
        this->hapcalls[1][seqnames[i]] = variantCalls();
        this->calls[seqnames[i]] = variantCalls();
    }

    // struct for storing each record
    rec = bcf_init();
    if (rec == NULL) {
        ERROR("failed to read VCF records");
        goto error2;
    }
    
    while (bcf_read(vcf, hdr, rec) == 0) {

        // new contig!
        seq = seqnames[rec->rid];
        if (rec->rid != prev_rid) {
            prev_rid = rec->rid;
            if (prev_rids.find(rec->rid) != prev_rids.end()) {
                ERROR("unsorted VCF, contig %s already parsed", seq.data());
            } else {
                INFO("parsing contig %s", seq.data());
                prev_end = {-g.gap*2, -g.gap*2};
                hap_var_idx = {0, 0};
                var_idx = 0;
            }
        }

        // unpack info (populates rec->d allele info)
        bcf_unpack(rec, BCF_UN_ALL);
        n++;

        // check that variant passed all filters
        pass = false;
        for (int i = 0; i < rec->d.n_flt; i++) {
            if (rec->d.flt[i] == pass_filter_id) pass = true;
        }
        if (!pass) continue;

        // parse GQ in either INT or FLOAT format, and GT
        if (int_qual) {
            ngq = bcf_get_format_int32(hdr, rec, "GQ", &gq, &ngq_arr);
            if (ngq == -2) {
                ngq = bcf_get_format_float(hdr, rec, "GQ", &fgq, &ngq_arr);
                gq[0] = int(fgq[0]);
                int_qual = false;
            }
        }
        else {
            ngq = bcf_get_format_float(hdr, rec, "GQ", &fgq, &ngq_arr);
            gq[0] = int(fgq[0]);
        }
        if ( ngq == -3 ) {
            if (!gq_warn) {
                WARN("no GQ tag at %s:%lu", seq.data(), rec->pos);
                gq_warn = true; // only warn once
            }
            gq[0] = 0;
        }
        ngt = bcf_get_format_int32(hdr, rec, "GT", &gt, &ngt_arr);
        if (ngt < 0) {
            ERROR("failed to read GT at %s:%lu\n", seq.data(), rec->pos);
        }

        // parse genotype info
        // gt[i]  gt[i] >> 1  gt[i] << 1  GT MEANING
        // 0      0           0           .  missing
        // 2      1           4           0  reference
        // 4      2           8           1  alternate1
        // 6      3           12          2  alternate2
        switch ( (gt[0] >> 1) + (gt[1] << 1) ) {
            case 0:  ngts[GT_DOT_DOT]++;   break;
            case 5:  ngts[GT_REF_REF]++;   break;
            case 6:  ngts[GT_ALT1_REF]++;  break;
            case 7:  ngts[GT_ALT2_REF]++;  break;
            case 9:  ngts[GT_REF_ALT1]++;  break;
            case 10: ngts[GT_ALT1_ALT1]++; break;
            case 11: ngts[GT_ALT2_ALT1]++; break;
            case 13: ngts[GT_REF_ALT2]++;  break;
            case 14: ngts[GT_ALT1_ALT2]++; break;
            case 15: ngts[GT_ALT2_ALT2]++; break;
            default: ngts[GT_OTHER]++;     break; // 3|1 etc
        }

        // parse variant type
        for (int hap = 0; hap < 2; hap++) {

            // get ref and allele, skipping ref calls
            std::string ref = rec->d.allele[0];
            int alt_idx = bcf_gt_allele(gt[hap]);
            if (alt_idx < 0) alt_idx = hap; // set 0|1 if .|.
            if (alt_idx == 0) continue; // nothing to do if reference
            std::string alt = rec->d.allele[alt_idx];

            // skip spanning deletion
            if (alt == "*") { ntypes[hap][TYPE_REF]++; continue; }

            // determine variant type
            int pos = rec->pos;
            int type = -1;
            int lm = 0; // match from left->right (trim prefix)
            int rm = -1;// match from right->left (simplify complex variants GRP->INDEL)
            int reflen = int(ref.size());
            int altlen = int(alt.size());
            if (altlen-reflen > 0) { // insertion
                while (lm < reflen && ref[lm] == alt[lm]) lm++;
                while (reflen+rm >= lm && 
                        ref[reflen+rm] == alt[altlen+rm]) rm--;
                if (lm > reflen+rm) type = TYPE_INS; else type = TYPE_GRP;
                pos += lm;
                alt = alt.substr(lm, altlen+rm-lm+2);
                ref = ref.substr(lm, reflen+rm-lm+2);

            } else if (altlen-reflen < 0) { // deletion
                while (lm < altlen && ref[lm] == alt[lm]) lm++;
                while (altlen+rm >= lm && 
                        ref[reflen+rm] == alt[altlen+rm]) rm--;
                if (lm > altlen+rm) type = TYPE_DEL; else type = TYPE_GRP;
                pos += lm;
                alt = alt.substr(lm, altlen+rm-lm+2);
                ref = ref.substr(lm, reflen+rm-lm+2);

            } else { // substitution
                if (ref.size() == 1) {
                    type = (ref[0] == alt[0] ? TYPE_REF : TYPE_SUB);
                    if (type == TYPE_REF) continue;
                }
                else {
                    if (ref.substr(1) == alt.substr(1)){
                        type = TYPE_SUB;
                        ref = ref[0]; alt = alt[0]; // chop off matches
                    }
                    else type = TYPE_GRP;
                }
            }

            // TODO: keep overlaps, test all non-overlapping subsets
            // skip overlapping variants
            if (prev_end[hap] > pos) { // warn if overlap
                WARN("potential overlap at %s:%i, skipping", seq.data(), pos);
                continue;
            }

            // calculate reference length of variant
            int rlen;
            switch (type) {
                case TYPE_INS:
                    rlen = 0; break;
                case TYPE_SUB:
                case TYPE_REF:
                    rlen = 1; break;
                case TYPE_DEL:
                case TYPE_GRP:
                    rlen = ref.size(); break;
                default:
                    ERROR("unexpected variant type (%d)", type);
                    break;
            }

            // add to all calls info
            if (pos - std::max(prev_end[0], prev_end[1]) > g.gap) {
                this->calls[seq].add_gap(var_idx);
                offs = {0, 0};
            }
            if (type == TYPE_GRP) { // split GRP into INS+DEL
                this->calls[seq].add_var(pos, 0, offs[0], offs[1], hap, TYPE_INS,
                        "", alt, gq[0], rec->qual); // add INS
                offs[hap] += altlen;
                this->calls[seq].add_var(pos, rlen, offs[0], offs[1], hap, TYPE_DEL,
                        ref, "", gq[0], rec->qual); // add DEL
                offs[hap] -= reflen;
                var_idx += 2;
            } else {
                this->calls[seq].add_var(pos, rlen, offs[0], offs[1], hap, type,
                        ref, alt, gq[0], rec->qual);
                offs[hap] += altlen - reflen;
                var_idx++;
            }

            // add to haplotype-specific calls info
            if (pos - prev_end[hap] > g.gap) {
                this->hapcalls[hap][seq].add_gap(hap_var_idx[hap]);
                hap_offs[hap] = 0;
            }
            if (type == TYPE_GRP) { // split GRP into INS+DEL
                this->hapcalls[hap][seq].add_var(pos, 0, hap_offs[0], // INS
                    hap_offs[1], hap, TYPE_INS, "", alt, gq[0], rec->qual);
                hap_offs[hap] += altlen;
                this->hapcalls[hap][seq].add_var(pos, rlen, hap_offs[0], // DEL
                    hap_offs[1], hap, TYPE_DEL, ref, "", gq[0], rec->qual);
                hap_offs[hap] -= reflen;
                hap_var_idx[hap] += 2;
            } else {
                this->hapcalls[hap][seq].add_var(pos, rlen, hap_offs[0], 
                        hap_offs[1], hap, type, ref, alt, gq[0], rec->qual);
                hap_offs[hap] += altlen - reflen;
                hap_var_idx[hap]++;
            }

            prev_end[hap] = pos + rlen;
            npass++;
            ntypes[hap][type]++;
        }
    }
    this->calls[seq].add_gap(var_idx);
    this->hapcalls[0][seq].add_gap(hap_var_idx[0]);
    this->hapcalls[1][seq].add_gap(hap_var_idx[1]);

    INFO("VCF contains %i sample(s) and %i records", bcf_hdr_nsamples(hdr), n);
    INFO("After splitting by haplotype, %i PASS all filters", npass);

    printf("\nSequence names:");
    for (int i = 0; i < nseq; i++) {
        if (i % 5 == 0) printf("\n  ");
        printf("[%2i] %s \t", i, seqnames[i]);
    }
    printf("\n");

    printf("\nGenotypes:\n");
    for (size_t i = 0; i < gt_strs.size(); i++) {
        printf("  %s  %i\n", gt_strs[i].data(), ngts[i]);
    }

    printf("\nVariant Types:\n");
    for (int h = 0; h < 2; h++) {
        printf("  Haplotype %i\n", h+1);
        for (size_t i = 0; i < type_strs.size(); i++) {
            printf("    %s  %i\n", type_strs[i].data(), ntypes[h][i]);
        }
    }

    if (g.print_verbosity >= 1) {
        printf("\nGroups:\n");
        for(int i = 0; i < nseq; i++) {
            if (this->calls[seqnames[i]].poss.size())  {
                printf("  Contig %s: %lu variants, %lu groups total\n", 
                        seqnames[i],
                        this->calls[seqnames[i]].poss.size(),
                        this->calls[seqnames[i]].gaps.size());
                for (int h = 0; h < 2; h++) {
                    printf("    Haplotype %i: %lu variants, %lu groups total\n", h+1,
                            this->hapcalls[h][seqnames[i]].poss.size(),
                            this->hapcalls[h][seqnames[i]].gaps.size());
                }
            }
        }
    }

    free(gq);
    free(fgq);
    free(gt);
    free(seqnames);
    bcf_hdr_destroy(hdr);
    bcf_close(vcf);
    bcf_destroy(rec);
    return;
error2:
    free(seqnames);
error1:
    bcf_close(vcf);
    bcf_hdr_destroy(hdr);
    return;

}
