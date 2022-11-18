#include <filesystem>

#include "htslib/vcf.h"

#include "globals.h"
#include "print.h"

void Globals::parse_args(int argc, char ** argv) {

    /* print help and exit */
    if (argc < 4) {
        if (argc == 2 && (std::string(argv[1]) == "-v" || 
                    std::string(argv[1]) == "--version")) {
            this->print_version();
            std::exit(0);
        }
        this->print_usage();
        std::exit(0);
    }

    if (argv[1][0] == '-' || argv[2][0] == '-' || argv[3][0] == '-') {
        WARN("Optional arguments should be provided AFTER mandatory arguments");
        this->print_usage();
    }

    /* verify input VCF/FASTA filepaths */
    this->query_vcf_fn = std::string(argv[1]);
    this->query_vcf_path = std::filesystem::path(this->query_vcf_fn);
    htsFile* query_vcf_fp = bcf_open(query_vcf_fn.data(), "r");
    if (query_vcf_fp == NULL) {
        ERROR("Failed to open query VCF file '%s'", query_vcf_fn.data());
    } else {
        bcf_close(query_vcf_fp);
    }

    this->truth_vcf_fn = std::string(argv[2]);
    this->truth_vcf_path = std::filesystem::path(this->truth_vcf_fn);
    htsFile* truth_vcf_fp = bcf_open(truth_vcf_fn.data(), "r");
    if (truth_vcf_fp == NULL) {
        ERROR("Failed to open truth VCF file '%s'", truth_vcf_fn.data());
    } else {
        bcf_close(truth_vcf_fp);
    }

    // load reference FASTA
    this->ref_fasta_fn = std::string(argv[3]);
    INFO("0. Loading reference FASTA '%s'", ref_fasta_fn.data());
    this->ref_fasta_fp = fopen(ref_fasta_fn.data(), "r");
    if (ref_fasta_fp == NULL) {
        ERROR("Failed to open reference FASTA file '%s'", ref_fasta_fn.data());
    }

    /* handle optional arguments */
    for (int i = 4; i < argc;) {
        if (std::string(argv[i]) == "-b" || std::string(argv[i]) == "--bed") {
            i++;
            if (i == argc) {
                ERROR("Option '-b' used without providing BED filename");
            }
            try {
                this->bed_fn = std::string(argv[i]);
                this->bed = bedData(std::string(argv[i++]));
                this->bed.check();
            } catch (const std::exception & e) {
                ERROR("Invalid BED filename provided");
            }
/*******************************************************************************/
        } else if (std::string(argv[i]) == "-r" || 
                std::string(argv[i]) == "--results-prefix") {
            i++;
            if (i == argc) {
                ERROR("Option '-r' used without providing prefix for storing results");
            }
            try {
                this->out_prefix = std::string(argv[i++]);
                std::filesystem::create_directory(
                        std::filesystem::path(out_prefix).parent_path());
            } catch (const std::exception & e) {
                ERROR("%s", e.what());
            }
/*******************************************************************************/
        } else if (std::string(argv[i]) == "-q" || 
                std::string(argv[i]) == "--min-qual") {
            i++;
            if (i == argc) {
                ERROR("Option '-q' used without providing minimum variant quality");
            }
            try {
                this->min_qual = std::stoi(argv[i++]);
            } catch (const std::exception & e) {
                ERROR("Invalid minimum variant quality provided");
            }
            if (g.min_qual < 0) {
                ERROR("Must provide non-negative minimum variant quality");
            }
/*******************************************************************************/
        } else if (std::string(argv[i]) == "-m" || 
                std::string(argv[i]) == "--max-qual") {
            i++;
            if (i == argc) {
                ERROR("Option '-m' used without providing maximum variant quality");
            }
            try {
                this->max_qual = std::stoi(argv[i++]);
            } catch (const std::exception & e) {
                ERROR("Invalid maximum variant quality provided");
            }
            if (g.max_qual < g.min_qual) {
                ERROR("Maximum variant quality must exceed minimum variant quality");
            }
/*******************************************************************************/
        } else if (std::string(argv[i]) == "-p" || 
                std::string(argv[i]) == "--print-verbosity") {
            i++;
            if (i == argc) {
                ERROR("Option '-p' used without providing print verbosity");
            }
            try {
                this->print_verbosity = std::stoi(argv[i++]);
            } catch (const std::exception & e) {
                ERROR("Invalid print verbosity provided");
            }
            if (this->print_verbosity < 0 || this->print_verbosity > 2) {
                ERROR("Print verbosity %d not a valid option (0,1,2)", 
                        this->print_verbosity);
            }
/*******************************************************************************/
        } else if (std::string(argv[i]) == "-h" || 
                std::string(argv[i]) == "--help") {
            i++;
            this->print_usage();
/*******************************************************************************/
        } else if (std::string(argv[i]) == "-v" || 
                std::string(argv[i]) == "--version") {
            i++;
            this->print_version();
/*******************************************************************************/
        } else if (std::string(argv[i]) == "-g" || 
                std::string(argv[i]) == "--supercluster-gap") {
            i++;
            if (i == argc) {
                ERROR("Option '-g' used without providing minimum gap between (super)clusters");
            }
            try {
                this->cluster_min_gap = std::stoi(argv[i++]);
            } catch (const std::exception & e) {
                ERROR("Invalid minimum gap between (super)clusters provided");
            }
            if (g.cluster_min_gap <= 0) {
                ERROR("Must provide positive minimum gap between (super)clusters");
            }

/*******************************************************************************/
        } else if (std::string(argv[i]) == "-s" || 
                std::string(argv[i]) == "--sub-penalty") {
            i++;
            if (i == argc) {
                ERROR("Option '-s' used without providing substitution penalty");
            } else if (this->query_penalties_set[PEN_SUB]) {
                ERROR("Query substitution penalty already set, cannot use '-s'");
            } else if (this->truth_penalties_set[PEN_SUB]) {
                ERROR("Truth substitution penalty already set, cannot use '-s'");
            }
            try {
                this->query_sub = std::stoi(argv[i]);
                this->query_penalties_set[PEN_SUB] = true;
                this->truth_sub = std::stoi(argv[i++]);
                this->truth_penalties_set[PEN_SUB] = true;
            } catch (const std::exception & e) {
                ERROR("Invalid substitution penalty provided");
            }
            if (this->query_sub < 0) {
                ERROR("Must provide non-negative substitution penalty");
            }
/*******************************************************************************/
        } else if (std::string(argv[i]) == "-cs" || 
                std::string(argv[i]) == "--query-sub-penalty") {
            i++;
            if (i == argc) {
                ERROR("Option '-cs' used without providing substitution penalty");
            } else if (this->query_penalties_set[PEN_SUB]) {
                ERROR("Query substitution penalty already set, cannot use '-cs'");
            }
            try {
                this->query_sub = std::stoi(argv[i++]);
                this->query_penalties_set[PEN_SUB] = true;
            } catch (const std::exception & e) {
                ERROR("Invalid query substitution penalty provided");
            }
            if (this->query_sub < 0) {
                ERROR("Must provide non-negative query substitution penalty");
            }
/*******************************************************************************/
        } else if (std::string(argv[i]) == "-ts" || 
                std::string(argv[i]) == "--truth-sub-penalty") {
            i++;
            if (i == argc) {
                ERROR("Option '-ts' used without providing substitution penalty");
            } else if (this->truth_penalties_set[PEN_SUB]) {
                ERROR("Truth substitution penalty already set, cannot use '-ts'");
            }
            try {
                this->truth_sub = std::stoi(argv[i++]);
                this->truth_penalties_set[PEN_SUB] = true;
            } catch (const std::exception & e) {
                ERROR("Invalid truth substitution penalty provided");
            }
            if (this->truth_sub < 0) {
                ERROR("Must provide non-negative truth substitution penalty");
            }
/*******************************************************************************/
        } else if (std::string(argv[i]) == "-o" || 
                std::string(argv[i]) == "--gap-open-penalty") {
            i++;
            if (i == argc) {
                ERROR("Option '-o' used without providing gap-opening penalty");
            } else if (this->query_penalties_set[PEN_OPEN]) {
                ERROR("Query gap-opening penalty already set, cannot use '-o'");
            } else if (this->truth_penalties_set[PEN_OPEN]) {
                ERROR("Truth gap-opening penalty already set, cannot use '-o'");
            }
            try {
                this->query_open = std::stoi(argv[i]);
                this->query_penalties_set[PEN_OPEN] = true;
                this->truth_open = std::stoi(argv[i++]);
                this->truth_penalties_set[PEN_OPEN] = true;
            } catch (const std::exception & e) {
                ERROR("Invalid gap-opening penalty provided");
            }
            if (this->query_open < 0) {
                ERROR("Must provide non-negative gap-opening penalty");
            }
/*******************************************************************************/
        } else if (std::string(argv[i]) == "-co" || 
                std::string(argv[i]) == "--query-gap-open-penalty") {
            i++;
            if (i == argc) {
                ERROR("Option '-co' used without providing gap-opening penalty");
            } else if (this->query_penalties_set[PEN_OPEN]) {
                ERROR("Query gap-opening penalty already set, cannot use '-co'");
            }
            try {
                this->query_open = std::stoi(argv[i++]);
                this->query_penalties_set[PEN_OPEN] = true;
            } catch (const std::exception & e) {
                ERROR("Invalid gap-opening penalty provided");
            }
            if (this->query_open < 0) {
                ERROR("Must provide non-negative gap-opening penalty");
            }
/*******************************************************************************/
        } else if (std::string(argv[i]) == "-to" || 
                std::string(argv[i]) == "--truth-gap-open-penalty") {
            i++;
            if (i == argc) {
                ERROR("Option '-to' used without providing gap-opening penalty");
            } else if (this->truth_penalties_set[PEN_OPEN]) {
                ERROR("Truth gap-opening penalty already set, cannot use '-to'");
            }
            try {
                this->truth_open = std::stoi(argv[i++]);
                this->truth_penalties_set[PEN_OPEN] = true;
            } catch (const std::exception & e) {
                ERROR("Invalid truth gap-opening penalty provided");
            }
            if (this->truth_open < 0) {
                ERROR("Must provide non-negative truth gap-opening penalty");
            }
/*******************************************************************************/
        } else if (std::string(argv[i]) == "-e" || 
                std::string(argv[i]) == "--gap-extend-penalty") {
            i++;
            if (i == argc) {
                ERROR("Option '-e' used without providing gap-extension penalty");
            } else if (this->query_penalties_set[PEN_EXTEND]) {
                ERROR("Query gap-extension penalty already set, cannot use '-e'");
            } else if (this->truth_penalties_set[PEN_EXTEND]) {
                ERROR("Truth gap-extension penalty already set, cannot use '-e'");
            }
            try {
                this->query_extend = std::stoi(argv[i]);
                this->query_penalties_set[PEN_EXTEND] = true;
                this->truth_extend = std::stoi(argv[i++]);
                this->truth_penalties_set[PEN_EXTEND] = true;
            } catch (const std::exception & e) {
                ERROR("Invalid gap-extension penalty provided");
            }
            if (this->query_extend < 0) {
                ERROR("Must provide non-negative gap-extension penalty");
            }
/*******************************************************************************/
        } else if (std::string(argv[i]) == "-ce" || 
                std::string(argv[i]) == "--query-gap-extend-penalty") {
            i++;
            if (i == argc) {
                ERROR("Option '-ce' used without providing gap-extension penalty");
            } else if (this->query_penalties_set[PEN_EXTEND]) {
                ERROR("Query gap-extension penalty already set, cannot use '-ce'");
            }
            try {
                this->query_extend = std::stoi(argv[i++]);
                this->query_penalties_set[PEN_EXTEND] = true;
            } catch (const std::exception & e) {
                ERROR("Invalid query gap-extension penalty provided");
            }
            if (this->query_extend < 0) {
                ERROR("Must provide non-negative query gap-extension penalty");
            }
/*******************************************************************************/
        } else if (std::string(argv[i]) == "-te" || 
                std::string(argv[i]) == "--truth-gap-extend-penalty") {
            i++;
            if (i == argc) {
                ERROR("Option '-te' used without providing gap-extension penalty");
            } else if (this->truth_penalties_set[PEN_EXTEND]) {
                ERROR("Truth gap-extension penalty already set, cannot use '-te'");
            }
            try {
                this->truth_extend = std::stoi(argv[i++]);
                this->truth_penalties_set[PEN_EXTEND] = true;
            } catch (const std::exception & e) {
                ERROR("Invalid truth gap-extension penalty provided");
            }
            if (this->truth_extend < 0) {
                ERROR("Must provide non-negative truth gap-extension penalty");
            }
/*******************************************************************************/
        } else if (std::string(argv[i]) == "-x" || 
                std::string(argv[i]) == "--exit-after-realign") {
            i++;
            g.exit = true;
/*******************************************************************************/
        } else if (std::string(argv[i]) == "--keep-truth") {
            i++;
            g.keep_truth = true;
/*******************************************************************************/
        } else if (std::string(argv[i]) == "--keep-query") {
            i++;
            g.keep_query = true;
/*******************************************************************************/
        } else if (std::string(argv[i]) == "--simple-cluster") {
            i++;
            g.simple_cluster = true;
        } else {
            ERROR("Unexpected option '%s'", argv[i]);
        }
    }
}

/* --------------------------------------------------------------------------- */

void Globals::print_version() const
{
    printf("%s v%s\n", this->PROGRAM.data(), this->VERSION.data());
}

void Globals::print_usage() const
{
    printf("Usage: vcfdist <query.vcf> <truth.vcf> <ref.fasta> [options]\n"); 

    printf("\nRequired:\n");
    printf("  <FILENAME>\tquery.vcf\tphased VCF file containing variant calls to evaluate \n");
    printf("  <FILENAME>\ttruth.vcf\tphased VCF file containing ground truth variant calls \n");
    printf("  <FILENAME>\tref.fasta\tFASTA file containing draft reference sequence \n");

    printf("\nOptions:\n");
    printf("  -b, --bed <FILENAME>\n");
    printf("      BED file containing regions to evaluate\n\n");

    printf("  -r, --results-prefix <FILENAME_PREFIX>\n");
    printf("      output filepath prefix (directories should contain trailing slashes)\n\n");

    printf("  -q, --min-qual <VALUE> [0]\n");
    printf("      minimum variant quality to be considered (lower qualities ignored)\n\n");

    printf("  -m, --max-qual <VALUE> [60]\n");
    printf("      maximum variant quality (higher qualities kept, thresholded)\n\n");

    printf("  -g, --supercluster-gap <VALUE> [50]\n");
    printf("      minimum base gap between independent superclusters\n\n");

    printf("  -s, --sub-penalty <VALUE> [1]\n");
    printf("      integer Smith-Waterman substitution penalty\n\n");

    printf("  -o, --gap-open-penalty <VALUE> [1]\n");
    printf("      integer Smith-Waterman gap opening penalty\n\n");

    printf("  -e, --gap-extend-penalty <VALUE> [1]\n");
    printf("      integer Smith-Waterman gap extension penalty\n\n");

    printf("  -p, --print-verbosity <VALUE> [0]\n");
    printf("      printing verbosity (0: default, 1: verbose, 2:debugging)\n\n");

    printf("  -x, --exit-after-realign\n");
    printf("      realign truth/query VCFs with Smith-Waterman params, then exit\n\n");

    printf("  -h, --help\n");
    printf("      show this help message\n\n");

    printf("  -v, --version\n");
    printf("      show version number (%s v%s)\n\n", this->PROGRAM.data(), this->VERSION.data());

    printf("\nAdvanced Options: (not recommended, only for evaluation)\n");
    printf("  --keep-query\n");
    printf("      do not realign query.vcf variants using Smith-Waterman\n\n");

    printf("  --keep-truth\n");
    printf("      do not realign truth.vcf variants using Smith-Waterman\n\n");

    printf("  --simple-cluster\n");
    printf("      instead of Smith-Waterman clustering, use gap-based clustering \n\n");

    printf("  -cs, --query-sub-penalty <VALUE> [1]\n");
    printf("      integer Smith-Waterman substitution penalty for query variants\n\n");

    printf("  -co, --query-gap-open-penalty <VALUE> [1]\n");
    printf("      integer Smith-Waterman gap opening penalty for query variants\n\n");

    printf("  -ce, --query-gap-extend-penalty <VALUE> [1]\n");
    printf("      integer Smith-Waterman gap extension penalty for query variants\n\n");

    printf("  -ts, --truth-sub-penalty <VALUE> [1]\n");
    printf("      integer Smith-Waterman substitution penalty for truth variants\n\n");

    printf("  -to, --truth-gap-open-penalty <VALUE> [1]\n");
    printf("      integer Smith-Waterman gap opening penalty for truth variants\n\n");

    printf("  -te, --truth-gap-extend-penalty <VALUE> [1]\n");
    printf("      integer Smith-Waterman gap extension penalty for truth variants\n\n");
}
