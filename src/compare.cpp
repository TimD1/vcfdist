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
#define TYPE_MNP 5
std::vector< std::string > type_strs { "REF", "SUB", "INS", "DEL", "GRP", "MNP" };

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
    /* std::vector<int> lens; */
    /* std::vector<int> inss; */
    /* std::vector<int> dels; */
};

class vcfData {
public:
    vcfData(htsFile* vcf) {

        // counters
        int nseq   = 0;   // number of sequences
        std::vector< std::vector<int> > 
            ntypes(2, std::vector<int>(type_strs.size(), 0));
        int n      = 0;   // total number of records in file
        int npass  = 0;   // records which PASS all filters

        // data
        int prev_end = 0; // previous variant end position

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
            this->calls[seqnames[i]] = variantCalls();
            this->calls[seqnames[i]].gaps.push_back(0);
        }

        // struct for storing each record
        rec = bcf_init();
        if (rec == NULL) {
            ERROR("failed to read VCF records");
            goto error2;
        }
        
        while (bcf_read(vcf, hdr, rec) == 0) {

            // unpack info (populates rec->d allele info)
            bcf_unpack(rec, BCF_UN_ALL);
            n++;

            // check that variant passed all filters
            pass = false;
            for (int i = 0; i < rec->d.n_flt; i++) {
                if (rec->d.flt[i] == pass_filter_id) pass = true;
            }
            if (!pass) continue;
            npass++;

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
                int allele_idx = bcf_gt_allele(gt[hap]);
                if (allele_idx < 0) allele_idx = hap; // set 0|1 if .|.
                if (allele_idx == 0) continue; // nothing to do if reference
                std::string allele = rec->d.allele[allele_idx];
                // skip spanning deletion
                if (allele == "*") { ntypes[hap][TYPE_REF]++; continue; }

                // determine variant type
                int indel_len = allele.size() - ref.size();
                int type = -1;
                if (indel_len > 0) { // insertion
                    type = TYPE_INS;
                } else if (indel_len < 0) { // deletion
                    type = TYPE_DEL;
                } else { // substitution
                    if (ref.size() == 1) type = TYPE_SUB;
                    else {
                        if (std::string(ref.data()+1) == std::string(allele.data()+1)){
                            type = TYPE_SUB;
                            ref = ref[0]; allele = allele[0]; // chop off matches
                        }
                        else type = TYPE_MNP;
                    }
                }
                ntypes[hap][type]++;

                if (rec->pos - prev_end > 50)
                    this->calls[seqnames[rec->rid]].gaps.push_back(npass-1);
                prev_end = rec->pos + rec->rlen;
                this->calls[seqnames[rec->rid]].poss.push_back(rec->pos);
                this->calls[seqnames[rec->rid]].haps.push_back(hap);
                this->calls[seqnames[rec->rid]].types.push_back(type);
                this->calls[seqnames[rec->rid]].refs.push_back(ref);
                this->calls[seqnames[rec->rid]].alts.push_back(allele);
                this->calls[seqnames[rec->rid]].gt_quals.push_back(gq[0]);
                this->calls[seqnames[rec->rid]].var_quals.push_back(rec->qual);
            }
        }
        this->calls[seqnames[rec->rid]].gaps.push_back(npass);

        fprintf(stderr, "VCF contains %i sample(s) and %i records, "
            "of which %i PASS all filters.\n", bcf_hdr_nsamples(hdr), n, npass);

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
        fprintf(stderr, "%i groups total\n", 
                int(this->calls[seqnames[rec->rid]].gaps.size()));

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
private:
    std::unordered_map<std::string, variantCalls> calls;
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
 
int main(int argc, char **argv) {

    // parse and store command-line args
    htsFile *ref_fasta, *calls_vcf, *truth_vcf;
    bedData * bed;
    std::string out_dir;
    Arguments args;
    int exit = parse_args(argc, argv, 
            ref_fasta, calls_vcf, truth_vcf, bed, out_dir);
    if (exit) std::exit(EXIT_SUCCESS);

    fprintf(stderr, "CALLS:\n");
    vcfData * calls = new vcfData(calls_vcf);

    /* fprintf(stderr, "\n\nTRUTH:\n"); */
    /* vcfData * truth = new vcfData(truth_vcf); */

        
    delete bed;
    delete calls;
    /* delete truth; */

    return EXIT_SUCCESS;
}
