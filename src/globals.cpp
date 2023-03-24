#include "htslib/vcf.h"

#include "globals.h"
#include "print.h"

void Globals::parse_args(int argc, char ** argv) {

    /* if required arguments are not provided, you can only print help and exit */
    bool print_help = false;
    if (argc < 4) {
        int i = 1;
        while (i < argc) {
            if (std::string(argv[i]) == "-a" || 
                    std::string(argv[i]) == "--advanced") {
                i++;
                print_help = true;
                g.advanced = true;
            } else if (std::string(argv[i]) == "-h" || 
                    std::string(argv[i]) == "--help") {
                i++;
                print_help = true;
            } else if (std::string(argv[i]) == "--version") {
                i++;
                this->print_version();
            }
        }

        if (print_help || argc == 1) this->print_usage();
        std::exit(0);
    }

    if (argv[1][0] == '-' || argv[2][0] == '-' || argv[3][0] == '-') {
        WARN("Optional arguments should be provided AFTER mandatory arguments, cannot use STDIN");
        this->print_usage();
    }

    /* verify input VCF/FASTA filepaths */
    this->query_vcf_fn = std::string(argv[1]);
    htsFile* query_vcf_fp = bcf_open(query_vcf_fn.data(), "r");
    if (query_vcf_fp == NULL) {
        ERROR("Failed to open query VCF file '%s'", query_vcf_fn.data());
    } else {
        bcf_close(query_vcf_fp);
    }

    this->truth_vcf_fn = std::string(argv[2]);
    htsFile* truth_vcf_fp = bcf_open(truth_vcf_fn.data(), "r");
    if (truth_vcf_fp == NULL) {
        ERROR("Failed to open truth VCF file '%s'", truth_vcf_fn.data());
    } else {
        bcf_close(truth_vcf_fp);
    }

    // load reference FASTA
    this->ref_fasta_fn = std::string(argv[3]);
    INFO("Loading reference FASTA '%s'", ref_fasta_fn.data());
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
                this->bed_exists = true;
                this->bed.check();
            } catch (const std::exception & e) {
                ERROR("Invalid BED filename provided");
            }
/*******************************************************************************/
        } else if (std::string(argv[i]) == "-p" || 
                std::string(argv[i]) == "--prefix") {
            i++;
            if (i == argc) {
                ERROR("Option '-p' used without providing prefix for storing results");
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
        } else if (std::string(argv[i]) == "-v" || 
                std::string(argv[i]) == "--verbosity") {
            i++;
            if (i == argc) {
                ERROR("Option '-v' used without providing printing verbosity");
            }
            try {
                this->verbosity = std::stoi(argv[i++]);
            } catch (const std::exception & e) {
                ERROR("Invalid printing verbosity provided");
            }
            if (this->verbosity < 0 || this->verbosity > 3) {
                ERROR("Printing verbosity %d not a valid option (0,1,2,3)", 
                        this->verbosity);
            }
/*******************************************************************************/
        } else if (std::string(argv[i]) == "-h" || 
                std::string(argv[i]) == "--help") {
            i++;
            print_help = true;
/*******************************************************************************/
        } else if (std::string(argv[i]) == "--version") {
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
        } else if (std::string(argv[i]) == "-x" || 
                std::string(argv[i]) == "--mismatch-penalty") {
            i++;
            if (i == argc) {
                ERROR("Option '-x' used without providing mismatch penalty");
            } else if (this->query_penalties_set[PEN_SUB]) {
                ERROR("Query mismatch penalty already set, cannot use '-x'");
            } else if (this->truth_penalties_set[PEN_SUB]) {
                ERROR("Truth mismatch penalty already set, cannot use '-x'");
            } else if (this->eval_penalties_set[PEN_SUB]) {
                ERROR("Eval mismatch penalty already set, cannot use '-x'");
            }
            try {
                this->eval_sub = std::stoi(argv[i]);
                this->eval_penalties_set[PEN_SUB] = true;
                this->query_sub = std::stoi(argv[i]);
                this->query_penalties_set[PEN_SUB] = true;
                this->truth_sub = std::stoi(argv[i++]);
                this->truth_penalties_set[PEN_SUB] = true;
            } catch (const std::exception & e) {
                ERROR("Invalid mismatch penalty provided");
            }
            if (this->query_sub < 0) {
                ERROR("Must provide non-negative mismatch penalty");
            }
/*******************************************************************************/
        } else if (std::string(argv[i]) == "-qs" || 
                std::string(argv[i]) == "--query-sub-penalty") {
            i++;
            if (i == argc) {
                ERROR("Option '-qs' used without providing substitution penalty");
            } else if (this->query_penalties_set[PEN_SUB]) {
                ERROR("Query substitution penalty already set, cannot use '-qs'");
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
        } else if (std::string(argv[i]) == "-es" || 
                std::string(argv[i]) == "--eval-sub-penalty") {
            i++;
            if (i == argc) {
                ERROR("Option '-es' used without providing substitution penalty");
            } else if (this->eval_penalties_set[PEN_SUB]) {
                ERROR("Eval substitution penalty already set, cannot use '-es'");
            }
            try {
                this->eval_sub = std::stoi(argv[i++]);
                this->eval_penalties_set[PEN_SUB] = true;
            } catch (const std::exception & e) {
                ERROR("Invalid eval substitution penalty provided");
            }
            if (this->eval_sub < 0) {
                ERROR("Must provide non-negative eval substitution penalty");
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
            } else if (this->eval_penalties_set[PEN_OPEN]) {
                ERROR("Eval gap-opening penalty already set, cannot use '-o'");
            }
            try {
                this->eval_open = std::stoi(argv[i]);
                this->eval_penalties_set[PEN_OPEN] = true;
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
        } else if (std::string(argv[i]) == "-qo" || 
                std::string(argv[i]) == "--query-gap-open-penalty") {
            i++;
            if (i == argc) {
                ERROR("Option '-qo' used without providing gap-opening penalty");
            } else if (this->query_penalties_set[PEN_OPEN]) {
                ERROR("Query gap-opening penalty already set, cannot use '-qo'");
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
        } else if (std::string(argv[i]) == "-eo" || 
                std::string(argv[i]) == "--eval-gap-open-penalty") {
            i++;
            if (i == argc) {
                ERROR("Option '-eo' used without providing gap-opening penalty");
            } else if (this->eval_penalties_set[PEN_OPEN]) {
                ERROR("Eval gap-opening penalty already set, cannot use '-eo'");
            }
            try {
                this->eval_open = std::stoi(argv[i++]);
                this->eval_penalties_set[PEN_OPEN] = true;
            } catch (const std::exception & e) {
                ERROR("Invalid eval gap-opening penalty provided");
            }
            if (this->eval_open < 0) {
                ERROR("Must provide non-negative eval gap-opening penalty");
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
            } else if (this->eval_penalties_set[PEN_EXTEND]) {
                ERROR("Eval gap-extension penalty already set, cannot use '-e'");
            }
            try {
                this->eval_extend = std::stoi(argv[i]);
                this->eval_penalties_set[PEN_EXTEND] = true;
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
        } else if (std::string(argv[i]) == "-qe" || 
                std::string(argv[i]) == "--query-gap-extend-penalty") {
            i++;
            if (i == argc) {
                ERROR("Option '-qe' used without providing gap-extension penalty");
            } else if (this->query_penalties_set[PEN_EXTEND]) {
                ERROR("Query gap-extension penalty already set, cannot use '-qe'");
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
        } else if (std::string(argv[i]) == "-ee" || 
                std::string(argv[i]) == "--eval-gap-extend-penalty") {
            i++;
            if (i == argc) {
                ERROR("Option '-ee' used without providing gap-extension penalty");
            } else if (this->eval_penalties_set[PEN_EXTEND]) {
                ERROR("Eval gap-extension penalty already set, cannot use '-ee'");
            }
            try {
                this->eval_extend = std::stoi(argv[i++]);
                this->eval_penalties_set[PEN_EXTEND] = true;
            } catch (const std::exception & e) {
                ERROR("Invalid eval gap-extension penalty provided");
            }
            if (this->eval_extend < 0) {
                ERROR("Must provide non-negative eval gap-extension penalty");
            }
/*******************************************************************************/
        } else if (std::string(argv[i]) == "-i" || 
                std::string(argv[i]) == "--max-cluster-iterations") {
            i++;
            if (i == argc) {
                ERROR("Option '-i' used without providing max cluster iterations");
            }
            try {
                this->max_cluster_itrs = std::stoi(argv[i++]);
            } catch (const std::exception & e) {
                ERROR("Invalid max cluster iterations provided");
            }
            if (this->max_cluster_itrs < 1) {
                ERROR("Max cluster iterations must be positive");
            }
/*******************************************************************************/
        } else if (std::string(argv[i]) == "-r" || 
                std::string(argv[i]) == "--realign-only") {
            i++;
            g.exit = true;
/*******************************************************************************/
        } else if (std::string(argv[i]) == "-a" || 
                std::string(argv[i]) == "--advanced") {
            i++;
            g.advanced = true;
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
/*******************************************************************************/
        } else {
            ERROR("Unexpected option '%s'", argv[i]);
        }
    }
    if (print_help) 
        this->print_usage();
    else {
        std::string cmd = "";
        for (int i = 0; i < argc; i++)
            cmd += " " + std::string(argv[i]);
        INFO("Command: '%s'", cmd.data()+1);
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

    printf("  -p, --prefix <FILENAME_PREFIX> [./]\n");
    printf("      output filepath prefix (directories should contain trailing slashes)\n\n");

    printf("  -q, --min-qual <VALUE> [%d]\n", g.min_qual);
    printf("      minimum variant quality to be considered (lower qualities ignored)\n\n");

    printf("  -m, --max-qual <VALUE> [%d]\n", g.max_qual);
    printf("      maximum variant quality (higher qualities kept, thresholded)\n\n");

    printf("  -i, --max-cluster-iterations<VALUE> [%d]\n", g.max_cluster_itrs);
    printf("      maximum iterations for growing clusters\n\n");

    printf("  -g, --supercluster-gap <VALUE> [%d]\n", g.cluster_min_gap);
    printf("      minimum base gap between independent superclusters\n\n");

    printf("  -x, --mismatch-penalty <VALUE> [%d]\n", g.eval_sub);
    printf("      integer Smith-Waterman substitution penalty\n\n");

    printf("  -o, --gap-open-penalty <VALUE> [%d]\n", g.eval_open);
    printf("      integer Smith-Waterman gap opening penalty\n\n");

    printf("  -e, --gap-extend-penalty <VALUE> [%d]\n", g.eval_extend);
    printf("      integer Smith-Waterman gap extension penalty\n\n");

    printf("  -v, --verbosity <VALUE> [%d]\n", g.verbosity);
    printf("      printing verbosity (0: default, 1: verbose, 2:debug, 3:verbose debug)\n\n");

    printf("  -r, --realign-only\n");
    printf("      realign truth/query VCFs with Smith-Waterman params, then exit\n\n");

    printf("  -h, --help\n");
    printf("      show this help message\n\n");

    printf("  -a, --advanced\n");
    printf("      show advanced options\n\n");

    printf("  --version\n");
    printf("      show version number (%s v%s)\n\n", 
            this->PROGRAM.data(), this->VERSION.data());

    if (!this->advanced) return;
    printf("\nAdvanced Options: (not recommended, only for evaluation)\n");
    printf("  --keep-query\n");
    printf("      do not realign query.vcf variants using Smith-Waterman\n\n");

    printf("  --keep-truth\n");
    printf("      do not realign truth.vcf variants using Smith-Waterman\n\n");

    printf("  --simple-cluster\n");
    printf("      instead of Smith-Waterman clustering, use gap-based clustering \n\n");

    printf("  -qs, --query-sub-penalty <VALUE> [%d]\n", g.query_sub);
    printf("      integer Smith-Waterman substitution penalty for query variants\n\n");

    printf("  -qo, --query-gap-open-penalty <VALUE> [%d]\n", g.query_open);
    printf("      integer Smith-Waterman gap opening penalty for query variants\n\n");

    printf("  -qe, --query-gap-extend-penalty <VALUE> [%d]\n", g.query_extend);
    printf("      integer Smith-Waterman gap extension penalty for query variants\n\n");

    printf("  -ts, --truth-sub-penalty <VALUE> [%d]\n", g.truth_sub);
    printf("      integer Smith-Waterman substitution penalty for truth variants\n\n");

    printf("  -to, --truth-gap-open-penalty <VALUE> [%d]\n", g.truth_open);
    printf("      integer Smith-Waterman gap opening penalty for truth variants\n\n");

    printf("  -te, --truth-gap-extend-penalty <VALUE> [%d]\n", g.truth_extend);
    printf("      integer Smith-Waterman gap extension penalty for truth variants\n\n");

    printf("  -es, --eval-sub-penalty <VALUE> [%d]\n", g.eval_sub);
    printf("      integer Smith-Waterman substitution penalty for evaluating distance\n\n");

    printf("  -eo, --eval-gap-open-penalty <VALUE> [%d]\n", g.eval_open);
    printf("      integer Smith-Waterman gap opening penalty for evaluating distance\n\n");

    printf("  -ee, --eval-gap-extend-penalty <VALUE> [%d]\n", g.eval_extend);
    printf("      integer Smith-Waterman gap extension penalty for evaluating distance\n\n");

}
