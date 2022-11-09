#include <string>
#include <unordered_map>
#include <vector>

#include "htslib/vcf.h"

#include "variant.h"
#include "print.h"
#include "dist.h"

/******************************************************************************/

ctgVariants::ctgVariants() { 
    this->n = 0; 
}

void ctgVariants::add_cluster(int g) { this->clusters.push_back(g); }

void ctgVariants::add_var(int pos, int rlen, uint8_t hap, uint8_t type, uint8_t loc,
        std::string ref, std::string alt, uint8_t orig_gt, float gq, float vq) {
    this->poss.push_back(pos);
    this->rlens.push_back(rlen);
    this->haps.push_back(hap);
    this->types.push_back(type);
    this->locs.push_back(loc);
    this->refs.push_back(ref);
    this->alts.push_back(alt);
    this->orig_gts.push_back(orig_gt);
    this->gt_quals.push_back(gq);
    this->var_quals.push_back(vq);
    this->n++;

    this->errtypes.push_back(ERRTYPE_UN);
    this->credit.push_back(0);
}

/******************************************************************************/

void variantData::write(std::string out_vcf_fn) {

    // VCF header
    FILE* out_vcf = fopen(out_vcf_fn.data(), "w");
    const std::chrono::time_point now{std::chrono::system_clock::now()};
    time_t tt = std::chrono::system_clock::to_time_t(now);
    tm local_time = *localtime(&tt);
    fprintf(out_vcf, "##fileformat=VCFv4.2\n");
    fprintf(out_vcf, "##fileDate=%04d%02d%02d\n", local_time.tm_year + 1900, 
            local_time.tm_mon + 1, local_time.tm_mday);
    for (size_t i = 0; i < this->contigs.size(); i++)
        fprintf(out_vcf, "##contig=<ID=%s,length=%d>\n", 
                this->contigs[i].data(), this->lengths[i]);
    fprintf(out_vcf, "##FILTER=<ID=PASS,Description=\"All filters passed\">\n");
    fprintf(out_vcf, "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n");
    fprintf(out_vcf, "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t%s\n",
            this->sample.data());

    // write variants
    for (std::string ctg : this->contigs) {
        std::vector<size_t> ptrs = {0, 0};
        while (ptrs[HAP1] < this->ctg_variants[HAP1][ctg]->poss.size() ||
                ptrs[HAP2] < this->ctg_variants[HAP2][ctg]->poss.size()) {

            // get next positions, set flags for which haps
            int pos_hap1 = ptrs[HAP1] < this->ctg_variants[HAP1][ctg]->poss.size() ? 
                this->ctg_variants[HAP1][ctg]->poss[ptrs[HAP1]] : std::numeric_limits<int>::max();
            int pos_hap2 = ptrs[HAP2] < this->ctg_variants[HAP2][ctg]->poss.size() ? 
                this->ctg_variants[HAP2][ctg]->poss[ptrs[HAP2]] : std::numeric_limits<int>::max();

            // indels include previous base, adjust position
            if (this->ctg_variants[HAP1][ctg]->types[ptrs[HAP1]] == TYPE_INS || 
                    this->ctg_variants[HAP1][ctg]->types[ptrs[HAP1]] == TYPE_DEL) pos_hap1--;
            if (this->ctg_variants[HAP2][ctg]->types[ptrs[HAP2]] == TYPE_INS || 
                    this->ctg_variants[HAP2][ctg]->types[ptrs[HAP2]] == TYPE_DEL) pos_hap2--;
            int pos = std::min(pos_hap1, pos_hap2);
            bool hap1 = (pos_hap1 == pos);
            bool hap2 = (pos_hap2 == pos);

            // add variants to output VCF file
            if (hap1 && hap2) {
                if (this->ctg_variants[HAP1][ctg]->refs[ptrs[HAP1]] == 
                        this->ctg_variants[HAP2][ctg]->refs[ptrs[HAP2]] &&
                        this->ctg_variants[HAP1][ctg]->alts[ptrs[HAP1]] == 
                        this->ctg_variants[HAP2][ctg]->alts[ptrs[HAP2]]) {
                    
                    // homozygous variant (1|1)
                    print_variant(out_vcf, ctg, pos, 
                            this->ctg_variants[HAP1][ctg]->types[ptrs[HAP1]],
                            this->ctg_variants[HAP1][ctg]->refs[ptrs[HAP1]],
                            this->ctg_variants[HAP1][ctg]->alts[ptrs[HAP1]],
                            this->ctg_variants[HAP1][ctg]->var_quals[ptrs[HAP1]], "1|1");
                    
                } else {
                    // two separate phased variants (0|1 + 1|0)
                    print_variant(out_vcf, ctg, pos, 
                            this->ctg_variants[HAP1][ctg]->types[ptrs[HAP1]],
                            this->ctg_variants[HAP1][ctg]->refs[ptrs[HAP1]],
                            this->ctg_variants[HAP1][ctg]->alts[ptrs[HAP1]],
                            this->ctg_variants[HAP1][ctg]->var_quals[ptrs[HAP1]], "1|0");
                    print_variant(out_vcf, ctg, pos, 
                            this->ctg_variants[HAP2][ctg]->types[ptrs[HAP2]],
                            this->ctg_variants[HAP2][ctg]->refs[ptrs[HAP2]],
                            this->ctg_variants[HAP2][ctg]->alts[ptrs[HAP2]],
                            this->ctg_variants[HAP2][ctg]->var_quals[ptrs[HAP2]], "0|1");
                }

            } else if (hap1) { // 1|0
                print_variant(out_vcf, ctg, pos, 
                        this->ctg_variants[HAP1][ctg]->types[ptrs[HAP1]],
                        this->ctg_variants[HAP1][ctg]->refs[ptrs[HAP1]],
                        this->ctg_variants[HAP1][ctg]->alts[ptrs[HAP1]],
                        this->ctg_variants[HAP1][ctg]->var_quals[ptrs[HAP1]], "1|0");

            } else if (hap2) { // 0|1
                print_variant(out_vcf, ctg, pos, 
                        this->ctg_variants[HAP2][ctg]->types[ptrs[HAP2]],
                        this->ctg_variants[HAP2][ctg]->refs[ptrs[HAP2]],
                        this->ctg_variants[HAP2][ctg]->alts[ptrs[HAP2]],
                        this->ctg_variants[HAP2][ctg]->var_quals[ptrs[HAP2]], "0|1");
            }

            // update pointers
            if (hap1) ptrs[HAP1]++;
            if (hap2) ptrs[HAP2]++;
        }
    }
    fclose(out_vcf);
}


void variantData::print_variant(FILE* out_fp, std::string ctg, int pos, int type,
        std::string ref, std::string alt, float qual, std::string gt) {

    char ref_base;
    switch (type) {
    case TYPE_SUB:
        fprintf(out_fp, "%s\t%d\t.\t%s\t%s\t%f\tPASS\t.\tGT\t%s\n", ctg.data(),
            pos+1, ref.data(), alt.data(), qual, gt.data());
        break;
    case TYPE_INS:
    case TYPE_DEL:
        ref_base = this->ref->fasta.at(ctg)[pos];
        fprintf(out_fp, "%s\t%d\t.\t%s\t%s\t%f\tPASS\t.\tGT\t%s\n", ctg.data(), 
                pos+1, (ref_base + ref).data(), (ref_base + alt).data(), 
                qual, gt.data());
        break;
    default:
        ERROR("print_variant not implemented for type %d", type);
    }
}

void variantData::set_header(const std::unique_ptr<variantData> & vcf) {
    this->sample = vcf->sample;
    this->contigs = vcf->contigs;
    this->lengths = vcf->lengths;
    this->ref = vcf->ref;
    for (std::string ctg : this->contigs)
        for (int hap = 0; hap < 2; hap++)
            this->ctg_variants[hap][ctg] = 
                std::shared_ptr<ctgVariants>(new ctgVariants());
}


/* Copy variants from `cigar` string to `variantData`. */
void variantData::add_variants(
        const std::vector<int> & cigar, 
        int hap, int ref_pos,
        const std::string & ctg, 
        const std::string & calls, 
        const std::string & ref) {

    int calls_idx = 0;
    int ref_idx = 0;
    for (size_t cig_idx = 0; cig_idx < cigar.size(); ) {
        int indel_len = 0;
        switch (cigar[cig_idx]) {

            case PTR_DIAG: // no variant, update pointers
                cig_idx += 2;
                ref_idx++;
                calls_idx++;
                break;

            case PTR_SUB: // substitution
                cig_idx += 2;
                this->ctg_variants[hap][ctg]->add_var(ref_pos+ref_idx, 1, hap, 
                        TYPE_SUB, INSIDE, std::string(1,ref[ref_idx]), 
                        std::string(1,calls[calls_idx]), GT_REF_REF, 50, 50);
                ref_idx++;
                calls_idx++;
                break;

            case PTR_DEL: // deletion
                cig_idx++; indel_len++;

                // multi-base deletion
                while (cig_idx < cigar.size() && cigar[cig_idx] == PTR_DEL) {
                    cig_idx++; indel_len++;
                }
                this->ctg_variants[hap][ctg]->add_var(ref_pos+ref_idx,
                        indel_len, hap, TYPE_DEL, INSIDE,
                        ref.substr(ref_idx, indel_len),
                        "", GT_REF_REF, 50, 50);
                ref_idx += indel_len;
                break;

            case PTR_INS: // insertion
                cig_idx++; indel_len++;

                // multi-base insertion
                while (cig_idx < cigar.size() && cigar[cig_idx] == PTR_INS) {
                    cig_idx++; indel_len++;
                }
                this->ctg_variants[hap][ctg]->add_var(ref_pos+ref_idx,
                        0, hap, TYPE_INS, INSIDE, "", 
                        calls.substr(calls_idx, indel_len), GT_REF_REF, 50, 50);
                calls_idx += indel_len;
                break;
        }
    }
}

/******************************************************************************/

variantData::variantData() : ctg_variants(2) { ; }

variantData::variantData(std::string vcf_fn, std::shared_ptr<fastaData> reference) : ctg_variants(2) {

    // set reference fasta pointer
    this->ref = reference;
    this->filename = vcf_fn;

    INFO(" ");
    INFO("Parsing VCF '%s'", vcf_fn.data());
    htsFile* vcf = bcf_open(vcf_fn.data(), "r");

    // counters
    int nseq   = 0;                     // number of sequences
    std::vector< std::vector<int> > 
        ntypes(2, std::vector<int>(type_strs.size(), 0));
    int n      = 0;                     // total number of records in file
    std::vector<int> npass  = {0, 0};   // records PASSing all filters

    // data
    std::vector<int> prev_end = {-g.cluster_min_gap*2, -g.cluster_min_gap*2};
    std::unordered_map<int, bool> prev_rids;
    int prev_rid = -1;
    std::unordered_map<int, int> seqlens;
    std::string seq;
    std::vector<int> nregions(region_strs.size(), 0);
    std::vector<int> pass_min_qual = {0, 0};

    // quality data for each call
    int ngq_arr = 0;
    int ngq     = 0;
    int * gq    = (int*) malloc(sizeof(int));
    float * fgq = (float*) malloc(sizeof(float));
    bool int_qual = true;
    int total_overlaps = 0;

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
    for (int i = 0; i < hdr->nhrec; i++) {

        /* // DEBUG HEADER PRINT */
        /* printf("%s=%s\n", hdr->hrec[i]->key, hdr->hrec[i]->value); */
        /* for (int j = 0; j < hdr->hrec[i]->nkeys; j++) { */
        /*     printf("  %s=%s\n", hdr->hrec[i]->keys[j], hdr->hrec[i]->vals[j]); */
        /* } */

        // search all FILTER lines
        if (hdr->hrec[i]->type == BCF_HL_FLT) {

            // select PASS filter
            bool is_pass_filter = false;
            for (int j = 0; j < hdr->hrec[i]->nkeys; j++)
                if (hdr->hrec[i]->keys[j] == std::string("ID") && 
                        hdr->hrec[i]->vals[j] == std::string("PASS"))
                    is_pass_filter = true;

            // save PASS filter index to keep only passing reads
            if (is_pass_filter)
                for (int j = 0; j < hdr->hrec[i]->nkeys; j++)
                    if (hdr->hrec[i]->keys[j] == std::string("IDX")) {
                        pass_filter_id = std::stoi(hdr->hrec[i]->vals[j]);
                        pass_found = true;
                    }
        }

        // store contig lengths (for output VCF)
        else if (hdr->hrec[i]->type == BCF_HL_CTG) {
            int length = -1;
            int idx = -1;
            for (int j = 0; j < hdr->hrec[i]->nkeys; j++) {
                if (hdr->hrec[i]->keys[j] == std::string("IDX"))
                    idx = std::stoi(hdr->hrec[i]->vals[j]);
                else if (hdr->hrec[i]->keys[j] == std::string("length"))
                    length = std::stoi(hdr->hrec[i]->vals[j]);
            }
            if (length >= 0 && idx >= 0) seqlens[idx] = length;
            else ERROR("VCF header contig line didn't have 'IDX' and 'length'");
        }
    }
    if (!pass_found)
        ERROR("Failed to find PASS FILTER in VCF");
    if (bcf_hdr_nsamples(hdr) != 1) 
        ERROR("Expected 1 sample, found %d", bcf_hdr_nsamples(hdr));
    this->sample = hdr->samples[0];

    // report names of all the sequences in the VCF file
    const char **seqnames = NULL;
    seqnames = bcf_hdr_seqnames(hdr, &nseq);
    if (seqnames == NULL) {
        ERROR("Failed to read VCF header");
        goto error1;
    }
    for(int i = 0; i < nseq; i++) {
        this->ctg_variants[HAP1][seqnames[i]] = std::shared_ptr<ctgVariants>(new ctgVariants());
        this->ctg_variants[HAP2][seqnames[i]] = std::shared_ptr<ctgVariants>(new ctgVariants());
    }

    // struct for storing each record
    rec = bcf_init();
    if (rec == NULL) {
        ERROR("Failed to read VCF records");
        goto error2;
    }
    
    while (bcf_read(vcf, hdr, rec) == 0) {

        // new contig!
        seq = seqnames[rec->rid];
        if (rec->rid != prev_rid) {

            // start new contig
            prev_rid = rec->rid;
            if (prev_rids.find(rec->rid) != prev_rids.end()) {
                ERROR("Unsorted VCF, contig %s already parsed", seq.data());
            } else {
                this->contigs.push_back(seq);
                this->lengths.push_back(seqlens[rec->rid]);
                prev_end = {-g.cluster_min_gap*2, -g.cluster_min_gap*2};
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

        // check that variant exceeds min_qual
        pass = rec->qual >= g.min_qual;
        pass_min_qual[pass]++;
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
                WARN("No GQ tag at %s:%lu", seq.data(), rec->pos);
                gq_warn = true; // only warn once
            }
            gq[0] = 0;
        }
        ngt = bcf_get_format_int32(hdr, rec, "GT", &gt, &ngt_arr);
        if (ngt < 0) {
            ERROR("Failed to read GT at %s:%lu\n", seq.data(), rec->pos);
        }

        // parse genotype info
        // gt[i]  gt[i] >> 1  gt[i] << 1  GT MEANING
        // 0      0           0           .  missing
        // 2      1           4           0  reference
        // 4      2           8           1  alternate1
        // 6      3           12          2  alternate2
        int orig_gt;
        switch ( (gt[0] >> 1) + ((gt[1] >> 1) << 2) ) {
            case 0:  orig_gt = GT_DOT_DOT;   break;
            case 5:  orig_gt = GT_REF_REF;   break;
            case 6:  orig_gt = GT_ALT1_REF;  break;
            case 7:  orig_gt = GT_ALT2_REF;  break;
            case 9:  orig_gt = GT_REF_ALT1;  break;
            case 10: orig_gt = GT_ALT1_ALT1; break;
            case 11: orig_gt = GT_ALT2_ALT1; break;
            case 13: orig_gt = GT_REF_ALT2;  break;
            case 14: orig_gt = GT_ALT1_ALT2; break;
            case 15: orig_gt = GT_ALT2_ALT2; break;
            default: orig_gt = GT_OTHER;     break; // 3|1 etc
        }
        ngts[orig_gt]++;

        // parse variant type
        for (int hap = 0; hap < 2; hap++) {

            // set simplified GT (0|1, 1|0, or 1|1), (0|0 skipped)
            int simple_gt = hap ? GT_REF_ALT1 : GT_ALT1_REF; // 0|1 or 1|0 default
            if (orig_gt == GT_ALT1_ALT1 || orig_gt == GT_ALT2_ALT2) // 1|1 overwrite
                simple_gt = GT_ALT1_ALT1;

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
                alt = alt.substr(lm, altlen+rm-lm+1);
                ref = ref.substr(lm, reflen+rm-lm+1);

            } else if (altlen-reflen < 0) { // deletion
                while (lm < altlen && ref[lm] == alt[lm]) lm++;
                while (altlen+rm >= lm && 
                        ref[reflen+rm] == alt[altlen+rm]) rm--;
                if (lm > altlen+rm) type = TYPE_DEL; else type = TYPE_GRP;
                pos += lm;
                alt = alt.substr(lm, altlen+rm-lm+1);
                ref = ref.substr(lm, reflen+rm-lm+1);

            } else { // substitution
                if (ref.size() == 1) {
                    type = (ref[0] == alt[0] ? TYPE_REF : TYPE_SUB);
                    if (type == TYPE_REF) continue;
                } else {
                    if (ref.substr(1) == alt.substr(1)){
                        type = TYPE_SUB;
                        ref = ref[0]; alt = alt[0]; // chop off matches
                    }
                    else type = TYPE_GRP;
                }
            }

            // TODO: keep overlaps, test all non-overlapping subsets?
            // skip overlapping variants
            if (prev_end[hap] > pos) { // warn if overlap
                if (g.print_verbosity >= 1) {
                    WARN("VCF variant overlap at %s:%i, skipping", seq.data(), pos);
                }
                total_overlaps++;
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
                    ERROR("Unexpected variant type: %d", type);
                    break;
            }

            // check that variant is in region of interest
            uint8_t loc = g.bed.contains(seq, pos, pos + rlen);
            switch (loc) {
                case OUTSIDE: 
                case OFF_CTG:
                    nregions[loc]++;
                    continue; // discard variant
                case INSIDE: 
                case BORDER:
                    nregions[loc]++;
                    break;
                default:
                    ERROR("Unexpected BED region type: %d", loc);
                    break;
            }

            // add to haplotype-specific calls info
            if (type == TYPE_GRP) { // split GRP into INS+DEL
                this->ctg_variants[hap][seq]->add_var(pos, 0, // INS
                    hap, TYPE_INS, loc, "", alt, simple_gt, gq[0], rec->qual);
                this->ctg_variants[hap][seq]->add_var(pos, rlen, // DEL
                    hap, TYPE_DEL, loc, ref, "", simple_gt, gq[0], rec->qual);
            } else {
                this->ctg_variants[hap][seq]->add_var(pos, rlen,
                        hap, type, loc, ref, alt, simple_gt, gq[0], rec->qual);
            }

            prev_end[hap] = pos + rlen;
            npass[hap]++;
            ntypes[hap][type]++;
        }
    }

    if (total_overlaps)
        WARN("%d total VCF variant overlaps skipped", total_overlaps);

    INFO("Contigs:");
    for (size_t i = 0; i < this->contigs.size(); i++) {
        INFO("  [%2lu] %s", i, this->contigs[i].data());
    }
    INFO(" ");

    INFO("Genotypes:");
    for (size_t i = 0; i < gt_strs.size(); i++) {
        INFO("  %s  %i", gt_strs[i].data(), ngts[i]);
    }
    INFO(" ");

    INFO("Variant exceeds Min Qual (%d):", g.min_qual);
    INFO("  FAIL  %d", pass_min_qual[FAIL]);
    INFO("  PASS  %d", pass_min_qual[PASS]);
    INFO(" ");

    INFO("Variants in BED Regions:");
    for (size_t i = 0; i < region_strs.size(); i++) {
        INFO("  %s  %i", region_strs[i].data(), nregions[i]);
    }
    INFO(" ");

    INFO("Variant Types:");
    for (int h = 0; h < 2; h++) {
        INFO("  Haplotype %i", h+1);
        for (size_t i = 0; i < type_strs.size(); i++) {
            INFO("    %s  %i", type_strs[i].data(), ntypes[h][i]);
        }
    }
    INFO(" ");

    if (g.print_verbosity >= 1) {
        INFO("Clusters:");
        for(int i = 0; i < nseq; i++) {
            for (int h = 0; h < 2; h++) {
                if (this->ctg_variants[h][seqnames[i]]->poss.size() > 0) {
                    INFO("    '%s' hap %i: %lu variants", 
                        this->contigs[i].data(),
                        h+1,
                        this->ctg_variants[h][seqnames[i]]->poss.size());
                }
            }
        }
        INFO(" ");
    }

    INFO("Overview:");
    INFO("  VARIANTS  %d", n);
    INFO("  KEPT HAP1 %d", npass[HAP1]);
    INFO("  KEPT HAP2 %d", npass[HAP2]);
    INFO(" ");

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
