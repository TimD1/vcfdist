#include <cstdio>
#include <cstdlib>
#include <string>
#include <exception>
#include <filesystem>
#include <unordered_map>
#include <vector>
#include <fstream>
#include <iostream>
#include <sstream>
#include <sys/time.h>

#include "vcf.h"
#include "vcfutils.h"

std::string VERSION = "0.0.1";
std::string PROGRAM = "vcfdist::compare";
void print_version() {
    printf("%s v%s\n", PROGRAM.data(), VERSION.data());
}

struct Arguments {
    std::string calls, truth, ref, bed, out;
    int gap = 50;
} args;

/* --------------------------------------------------------------------------- */

#define WARN(f_, ...)                                           \
{                                                               \
    struct tm _tm123_;                                          \
    struct timeval _xxtv123_;                                   \
    gettimeofday(&_xxtv123_, NULL);                             \
    localtime_r(&_xxtv123_.tv_sec, &_tm123_);                   \
    fprintf(stderr, "\033[33m[WARN %s %02d:%02d:%02d]\033[0m ", \
            PROGRAM.data(), _tm123_.tm_hour,_tm123_.tm_min,     \
            _tm123_.tm_sec);                                    \
    fprintf(stderr, (f_), ##__VA_ARGS__);                       \
    fprintf(stderr, "\n");                                      \
};

#define INFO(f_, ...)                                           \
{                                                               \
    struct tm _tm123_;                                          \
    struct timeval _xxtv123_;                                   \
    gettimeofday(&_xxtv123_, NULL);                             \
    localtime_r(&_xxtv123_.tv_sec, &_tm123_);                   \
    fprintf(stderr, "\033[32m[INFO %s %02d:%02d:%02d]\033[0m ", \
            PROGRAM.data(), _tm123_.tm_hour,_tm123_.tm_min,     \
            _tm123_.tm_sec);                                    \
    fprintf(stderr, (f_), ##__VA_ARGS__);                       \
    fprintf(stderr, "\n");                                      \
};

#define ERROR(f_, ...)                                           \
{                                                                \
    struct tm _tm123_;                                           \
    struct timeval _xxtv123_;                                    \
    gettimeofday(&_xxtv123_, NULL);                              \
    localtime_r(&_xxtv123_.tv_sec, &_tm123_);                    \
    fprintf(stderr, "\033[31m[ERROR %s %02d:%02d:%02d]\033[0m ", \
            PROGRAM.data(), _tm123_.tm_hour,_tm123_.tm_min,      \
            _tm123_.tm_sec);                                     \
    fprintf(stderr, (f_), ##__VA_ARGS__);                        \
    fprintf(stderr, "\n");                                       \
};

/* --------------------------------------------------------------------------- */

void print_usage() {
    printf("Usage: vcfdist compare <calls.vcf> <truth.vcf> <ref.fasta> [options]\n"); 
    printf("Required:\n");
    printf("  FILE\tcalls.vcf\tVCF file containing variant calls to evaluate \n");
    printf("  FILE\ttruth.vcf\tVCF file containing ground truth variant calls \n");
    printf("  FILE\tref.fasta\tFASTA file containing reference sequence \n");
    printf("Options:\n");
    printf("  -b FILE\tBED file containing regions to evaluate\n");
    printf("  -o DIR\toutput directory\n");
    printf("  -g GAP size for independent groups [50]\n");
    printf("  -h\t\tshow this help message\n");
    printf("  -v\t\tshow version number\n");
}

/* --------------------------------------------------------------------------- */

struct contigRegions {
    std::vector<int> starts;
    std::vector<int> stops;
};

class bedData {
public:

    // constructors
    bedData() = default;
    bedData(std::string bed_fn) {
        std::ifstream bed(bed_fn);
        std::string region;
        while (getline(bed, region)) {
            std::stringstream ss(region);
            std::string contig, start, stop;
            getline(ss, contig, '\t');
            getline(ss, start, '\t');
            getline(ss, stop, '\t');
            this->add(contig, std::stoi(start), std::stoi(stop));
        }
    }

    void add(const std::string & contig, const int & start, const int & stop) {
        if (regions.find(contig) == regions.end()) {
            regions[contig] = contigRegions();
        }
        regions[contig].starts.push_back(start);
        regions[contig].stops.push_back(stop);
    }

    bool includes(std::string contig, long pos) {
        return false;
    }

    operator std::string() const {
        std::string bed_regions = "";
        for (const auto &[contig, region_list]: this->regions) {
            bed_regions += contig + ":\n";
            for (size_t i = 0; i < region_list.starts.size(); i++) {
                bed_regions += "\t" + std::to_string(region_list.starts[i]) + "-" + \
                               std::to_string(region_list.stops[i]) + "\n";
            }
        }
        return bed_regions;
    }

private:
    std::unordered_map<std::string, contigRegions> regions;
};

/* --------------------------------------------------------------------------- */

#define TYPE_REF 0
#define TYPE_SUB 1
#define TYPE_INS 2
#define TYPE_DEL 3
#define TYPE_GRP 4
std::vector< std::string > type_strs { "REF", "SUB", "INS", "DEL", "GRP"};

#define GT_DOT_DOT   0
#define GT_REF_REF   1
#define GT_REF_ALT1  2
#define GT_REF_ALT2  3
#define GT_ALT1_REF  4
#define GT_ALT1_ALT1 5
#define GT_ALT1_ALT2 6
#define GT_ALT2_REF  7
#define GT_ALT2_ALT1 8
#define GT_ALT2_ALT2 9
#define GT_OTHER     10
std::vector< std::string > gt_strs {
".|.", "0|0", "0|1", "0|2", "1|0", "1|1", "1|2", "2|0", "2|1", "2|2", "?|?"
};

struct variantCalls {
    std::vector<int> gaps;          // indices of gaps in this struct's vectors
    std::vector<int> poss;          // variant start positions
    std::vector<uint8_t> haps;      // variant haplotype
    std::vector<uint8_t> types;     // variant type: NONE, SUB, INS, DEL, GRP
    std::vector<std::string> refs;  // variant reference allele
    std::vector<std::string> alts;  // variant alternate allele (always one)
    std::vector<float> var_quals;   // variant quality (0-60)
    std::vector<float> gt_quals;    // genotype quality (0-60)
    std::vector<int> rlens;          // reference lengths
    /* std::vector<int> inss; */
    /* std::vector<int> dels; */
};

class vcfData {
public:
    vcfData(htsFile* vcf) : hapcalls(2) {

        // counters
        int nseq   = 0;                     // number of sequences
        std::vector< std::vector<int> > 
            ntypes(2, std::vector<int>(type_strs.size(), 0));
        int n      = 0;                     // total number of records in file
        int npass  = 0;                     // records PASSing all filters

        // data
        std::vector<int> prev_end = {-args.gap*2, -args.gap*2};
        std::unordered_map<int, bool> prev_rids;
        int prev_rid = -1;
        std::vector<int> var_idx  = {0, 0};   // indices per hap/contig

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
            if (rec->rid != prev_rid) {
                prev_rid = rec->rid;
                if (prev_rids.find(rec->rid) != prev_rids.end()) {
                    ERROR("unsorted VCF, contig %s already parsed", 
                            seqnames[rec->rid]);
                } else {
                    INFO("contig %s", seqnames[rec->rid]);
                    var_idx = {0, 0};
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
                    WARN("no GQ tag at %s:%i", seqnames[rec->rid], int(rec->pos));
                    gq_warn = true; // only warn once
                }
                gq[0] = 0;
            }
            ngt = bcf_get_format_int32(hdr, rec, "GT", &gt, &ngt_arr);
            if (ngt < 0) {
                ERROR("failed to read GT at %s:%i\n", 
                        seqnames[rec->rid], int(rec->pos));
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
                int indel_len = alt.size() - ref.size();
                int type = -1;
                int pos = rec->pos;
                size_t matches = 0;
                if (indel_len > 0) { // insertion
                    while (matches < ref.size() && ref[matches] == alt[matches]) matches++;
                    if (matches == ref.size()) type = TYPE_INS; else type = TYPE_GRP;
                    pos += matches;
                    alt = alt.substr(matches);
                    ref = ref.substr(matches);
                } else if (indel_len < 0) { // deletion
                    while (matches < alt.size() && ref[matches] == alt[matches]) matches++;
                    if (matches == alt.size()) type = TYPE_DEL; else type = TYPE_GRP;
                    pos += matches;
                    ref = ref.substr(matches);
                    alt = alt.substr(matches);
                } else { // substitution
                    if (ref.size() == 1) type = TYPE_SUB;
                    else {
                        if (ref.substr(1) == alt.substr(1)){
                            type = TYPE_SUB;
                            ref = ref[0]; alt = alt[0]; // chop off matches
                        }
                        else type = TYPE_GRP;
                    }
                }

                int rlen;
                switch (type) {
                    case TYPE_INS:
                        rlen = 0; break;
                    case TYPE_SUB:
                        rlen = 1; break;
                    case TYPE_DEL:
                    case TYPE_GRP:
                        rlen = ref.size(); break;
                    default:
                        ERROR("unexpected variant type (%d)", type);
                        break;
                }

                // add to all calls info
                if (pos - std::max(prev_end[0], prev_end[1]) > args.gap)
                    this->calls[seqnames[rec->rid]].gaps
                        .push_back(var_idx[0]+var_idx[1]);
                this->calls[seqnames[rec->rid]].poss.push_back(pos);
                this->calls[seqnames[rec->rid]].rlens.push_back(rlen);
                this->calls[seqnames[rec->rid]].haps.push_back(hap);
                this->calls[seqnames[rec->rid]].types.push_back(type);
                this->calls[seqnames[rec->rid]].refs.push_back(ref);
                this->calls[seqnames[rec->rid]].alts.push_back(alt);
                this->calls[seqnames[rec->rid]].gt_quals.push_back(gq[0]);
                this->calls[seqnames[rec->rid]].var_quals.push_back(rec->qual);

                // add to haplotype-specific calls info
                if (pos - prev_end[hap] > args.gap)
                    this->hapcalls[hap][seqnames[rec->rid]].gaps
                        .push_back(var_idx[hap]);
                this->hapcalls[hap][seqnames[rec->rid]].poss.push_back(pos);
                this->hapcalls[hap][seqnames[rec->rid]].rlens.push_back(rlen);
                this->hapcalls[hap][seqnames[rec->rid]].haps.push_back(hap);
                this->hapcalls[hap][seqnames[rec->rid]].types.push_back(type);
                this->hapcalls[hap][seqnames[rec->rid]].refs.push_back(ref);
                this->hapcalls[hap][seqnames[rec->rid]].alts.push_back(alt);
                this->hapcalls[hap][seqnames[rec->rid]].gt_quals.push_back(gq[0]);
                this->hapcalls[hap][seqnames[rec->rid]].var_quals.push_back(rec->qual);

                npass++;
                ntypes[hap][type]++;
                var_idx[hap]++;
                prev_end[hap] = pos + rlen;
            }
        }
        this->calls[seqnames[rec->rid]].gaps.push_back(var_idx[0]+var_idx[1]);
        this->hapcalls[0][seqnames[rec->rid]].gaps.push_back(var_idx[0]);
        this->hapcalls[1][seqnames[rec->rid]].gaps.push_back(var_idx[1]);

        fprintf(stderr, "VCF contains %i sample(s) and %i records, "
            "of which %i PASS all filters.\n", 
            bcf_hdr_nsamples(hdr), n, npass);

        fprintf(stderr, "\nSequence names:");
        for (int i = 0; i < nseq; i++) {
            if (i % 5 == 0) fprintf(stderr, "\n  ");
            fprintf(stderr, "[%2i] %s \t", i, seqnames[i]);
        }
        fprintf(stderr, "\n");

        fprintf(stderr, "\nGenotypes:\n");
        for (size_t i = 0; i < gt_strs.size(); i++) {
            fprintf(stderr, "  %s  %i\n", gt_strs[i].data(), ngts[i]);
        }

        fprintf(stderr, "\nVariant Types:\n");
        for (int h = 0; h < 2; h++) {
            fprintf(stderr, "  Haplotype %i\n", h+1);
            for (size_t i = 0; i < type_strs.size(); i++) {
                fprintf(stderr, "    %s  %i\n", type_strs[i].data(), ntypes[h][i]);
            }
        }

        fprintf(stderr, "\nGroups:\n");
        for(int i = 0; i < nseq; i++) {
            if (this->calls[seqnames[i]].poss.size())  {
                fprintf(stderr, "  Contig %s: %lu variants, %lu groups total\n", 
                        seqnames[i],
                        this->calls[seqnames[i]].poss.size(),
                        this->calls[seqnames[i]].gaps.size());
                for (int h = 0; h < 2; h++) {
                    fprintf(stderr, "    Haplotype %i: %lu variants, %lu groups total\n", h+1,
                            this->hapcalls[h][seqnames[i]].poss.size(),
                            this->hapcalls[h][seqnames[i]].gaps.size());
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
    std::unordered_map<std::string, variantCalls> calls;
    std::vector< std::unordered_map<std::string, variantCalls> > hapcalls;
};

/* --------------------------------------------------------------------------- */

int parse_args(
        int argc, 
        char ** argv, 
        htsFile *& ref_fasta, 
        htsFile *& calls_vcf, 
        htsFile *& truth_vcf, 
        bedData *& bed,
        std::string & out_dir
        ) {

    /* print help and exit */
    if (argc < 4) {
        if (argc == 2 && std::string(argv[1]) == "-v") {
            print_version();
            return 1;
        }
        print_usage();
        return 1;
    }

    /* verify input VCF/FASTA filepaths */
    std::string ref_fasta_fn, calls_vcf_fn, truth_vcf_fn;
    calls_vcf_fn = std::string(argv[1]);
    calls_vcf = bcf_open(calls_vcf_fn.data(), "r");
    if (calls_vcf == NULL) {
        ERROR("failed to open calls_vcf file '%s'.", calls_vcf_fn.data());
        return 1;
    }
    args.calls = calls_vcf_fn;
    truth_vcf_fn = std::string(argv[2]);
    truth_vcf = bcf_open(truth_vcf_fn.data(), "r");
    if (truth_vcf == NULL) {
        ERROR("failed to open truth_vcf file '%s'.", truth_vcf_fn.data());
        return 1;
    }
    args.truth = truth_vcf_fn;
    ref_fasta_fn = std::string(argv[3]);
    ref_fasta = bcf_open(ref_fasta_fn.data(), "r");
    if (ref_fasta == NULL) {
        ERROR("failed to open ref_fasta file '%s'.", ref_fasta_fn.data());
        return 1;
    }
    args.ref = ref_fasta_fn;

    /* handle optional arguments */
    for (int i = 4; i < argc;) {
        if (std::string(argv[i]) == "-b") {
            i++;
            if (i == argc) {
                ERROR("option '-b' used without providing BED.");
                return 1;
            }
            try {
                args.bed = std::string(argv[i]);
                bed = new bedData(std::string(argv[i++]));
            } catch (const std::exception & e) {
                ERROR("%s", e.what());
                return 1;
            }
        }
        else if (std::string(argv[i]) == "-o") {
            i++;
            if (i == argc) {
                ERROR("option '-o' used without providing DIR.");
                return 1;
            }
            try {
                args.out = argv[i];
                out_dir = std::string(argv[i++]);
                std::filesystem::create_directory(out_dir);
            } catch (const std::exception & e) {
                ERROR("%s", e.what());
                return 1;
            }
        }
        else if (std::string(argv[i]) == "-h") {
            i++;
            print_usage();
            return 1;
        }
        else if (std::string(argv[i]) == "-v") {
            i++;
            print_version();
            return 1;
        }
        else if (std::string(argv[i]) == "-g") {
            i++;
            args.gap = std::stoi(argv[i++]);
        }
        else {
            ERROR("unexpected option '%s'.", argv[i]);
            return 1;
        }
    }

    return 0;
}

/* --------------------------------------------------------------------------- */

int ed_align(vcfData* & vcf, htsFile* ref_fasta) {

    // iterate over each haplotype
    for (int h = 0; h < 2; h++) {

        // iterate over each contig
        for (auto itr = vcf->hapcalls[h].begin(); 
                itr != vcf->hapcalls[h].end(); itr++) {
            std::string ctg = itr->first;
            variantCalls vars = itr->second;
            if (vars.poss.size() == 0) continue;
            fprintf(stderr, "\nhap %i ctg %s: %lu variants\n", h, 
                    ctg.data(), vars.poss.size());

            // iterate over each group of variants
            for (size_t i = 0; i < vars.gaps.size()-1; i++) {
                int beg_idx = vars.gaps[i];
                int end_idx = vars.gaps[i+1];
                int beg = vars.poss[beg_idx]-1;
                int end = vars.poss[end_idx-1] + vars.rlens[end_idx-1]+1;
                int subs = 0;
                int inss = 0;
                int dels = 0;
                int grps = 0;
                const char* ref;
                const char* alt;

                // iterate over variants, summing edit distance
                for (int j = beg_idx; j < end_idx; j++) {
                    switch (vars.types[j]) {
                        case TYPE_SUB: subs++; break;
                        case TYPE_INS: inss += vars.alts[j].size(); break;
                        case TYPE_DEL: dels += vars.refs[j].size(); break;
                        case TYPE_GRP: grps += std::max(vars.refs[j].size(), 
                                               vars.alts[j].size()); break;
                        default: ERROR("unexpected type") std::exit(1); break;
                    }
                }

                /* if (grps > 0) { */
                if (end_idx-beg_idx > 1 && end-beg < 50) {
                    fprintf(stderr, "  group %i: %d variants, range %d\t(%d-%d),\tED %d (%dS %dI %dD %dG)\n",
                            int(i), end_idx-beg_idx, end-beg, beg, end,
                            subs*2 + inss + dels + grps, subs, inss, dels, grps);
                    for (int j = beg_idx; j < end_idx; j++) {
                        if (vars.refs[j].size() == 0) ref = "-"; else ref = vars.refs[j].data();
                        if (vars.alts[j].size() == 0) alt = "-"; else alt = vars.alts[j].data();
                        fprintf(stderr, "    %s:%i %s\t%s\t%s\n", ctg.data(), vars.poss[j],
                                type_strs[vars.types[j]].data(), ref, alt);
                    }
                }
            }
        }
    }

    return 0;

}

/* --------------------------------------------------------------------------- */
 
int main(int argc, char **argv) {

    // parse and store command-line args
    htsFile *ref_fasta, *calls_vcf, *truth_vcf;
    bedData * bed;
    std::string out_dir;
    Arguments args;
    int exit = parse_args(argc, argv, 
            ref_fasta, calls_vcf, truth_vcf, bed, out_dir);
    if (exit) std::exit(EXIT_SUCCESS);

    /* fprintf(stderr, "CALLS:\n"); */
    /* vcfData * calls = new vcfData(calls_vcf); */
    /* ed_align(calls, ref_fasta); */
    /* delete calls; */

    fprintf(stderr, "\n\nTRUTH:\n");
    vcfData * truth = new vcfData(truth_vcf);
    ed_align(truth, ref_fasta);
    delete truth;
        
    delete bed;

    return EXIT_SUCCESS;
}
