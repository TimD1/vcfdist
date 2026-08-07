/**
 * @file globals.cpp
 * @brief Global configuration object, string table definitions, and utility functions.
 */
#include "htslib/vcf.h"

#include <sstream>
#include <sys/stat.h>

#include "globals.h"
#include "print.h"
#include "timer.h"

Globals g;
/** @brief Program version string. */
const std::string Globals::VERSION = "3.0.0-b0";
/** @brief Program name string. */
const std::string Globals::PROGRAM = "vcfdist";
/** @brief String representations of QUERY/TRUTH callset indices. */
std::vector<std::string> callset_strs = {"QUERY", "TRUTH"};
/** @brief String representations of ERRTYPE_* constants (TP, FP, FN, unknown). */
std::vector<std::string> error_strs = {"TP", "FP", "FN", "??"};
/** @brief String representations of GT_* genotype constants. */
std::vector<std::string> gt_strs =
    {"0", "1", "0|0", "0|1", "1|0", "1|1", "1|2", "2|1", ".|.", "X|.", "X|Y" };
/** @brief String representations of PHASE_* constants (keep, swap, missing). */
std::vector<std::string> phase_strs = {"0", "1", "."};
/** @brief String representations of AC_ERR_* allele count error types. */
std::vector<std::string> ac_strs = {".", ".", ".", ".", "+", ".", "-", ".", "."};
/** @brief String representations of BED_* location constants. */
std::vector<std::string> region_strs = {"OUTSIDE", "INSIDE", "BORDER", "OFF_CTG"};
/** @brief String representations of SWITCHTYPE_* switch/flip error type constants. */
std::vector<std::string> switch_strs =
    {"FLIP", "SWITCH", "SWITCH+FLIP", "SWITCH_ERR", "FLIP_BEG", "FLIP_END", "NONE"};
/** @brief String names for pipeline stage timers in TIME_* index order. */
std::vector<std::string> timer_strs =
    {"reading", "clustering", "alignment eval", "phasing", "writing", "total"};
/** @brief String representations of TYPE_* variant type constants. */
std::vector<std::string> type_strs = {"REF", "SNP", "INS", "DEL", "CPX"};
/** @brief String representations of VARTYPE_* size-class constants. */
std::vector<std::string> vartype_strs = {"SNP", "INDEL", "SV", "ALL"};

/**
 * @brief Parses command-line arguments and initializes global configuration.
 * @param[in] argc Argument count
 * @param[in] argv Argument vector
 * @note Required args: query.vcf, truth.vcf, ref.fasta (must be first 3). Optional flag groups:
 *       input/output (-b, -v, -p, -n), variant filtering (-f, -l, -sv, -q, -mq),
 *       clustering (-s), precision-recall (-ct), resources (-t, -r), misc (-h, -ci).
 * @throws Errors on invalid file paths, out-of-range parameters, or format errors.
 */
void Globals::parse_args(int argc, char ** argv) {

    /* if required arguments are not provided, you can only print help and exit */
    bool print_cite = false;
    bool print_help = false;
    if (argc < 4) {
        int i = 1;
        while (i < argc) {
            if (std::string(argv[i]) == "-h" || 
                    std::string(argv[i]) == "--help") {
                i++;
                print_help = true;
            } else if (std::string(argv[i]) == "-v" ||
                    std::string(argv[i]) == "--version") {
                i++;
                this->print_version();
            } else if (std::string(argv[i]) == "-ci" || 
                    std::string(argv[i]) == "--citation") {
                i++;
                print_cite = true;
            } else {
                print_help = true;
                WARN("Invalid usage.");
                break;
            }
        }

        if (print_help || argc == 1) this->print_usage();
        if (print_cite) this->print_citation();
        std::exit(0);
    }

    // parse verbosity first
    for (int i = 0; i+1 < argc; i++) {
        if (std::string(argv[i]) == "-v" || 
                std::string(argv[i]) == "--verbosity") {
            i++;
            if (i == argc) {
                ERROR("Option '--verbosity' used without providing printing verbosity");
            }
            try {
                this->verbosity = std::stoi(argv[i]);
            } catch (const std::exception & e) {
                ERROR("Invalid printing verbosity provided");
            }
            if (this->verbosity < 0 || this->verbosity > 2) {
                ERROR("Printing verbosity %d not a valid option (0,1,2)", 
                        this->verbosity);
            }
            break;
        }
    }

    this->cmd = std::string(argv[0]);
    for (int i = 1; i < argc; i++)
        this->cmd += " " + std::string(argv[i]);
    if (this->verbosity >= 1) INFO("Command: '%s'", this->cmd.data());

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
    if (this->verbosity >= 1) {
        INFO(" ");
        INFO("%s[%d/%d] Loading reference FASTA%s '%s'", COLOR_PURPLE,
                int(idx(TIME_READ)), int(idx(TIME_TOTAL))-1, COLOR_WHITE, ref_fasta_fn.data());
    }
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
            } catch (const std::exception & e) {
                ERROR("%s", e.what());
            }
/**************************************************************************************************/
        } else if (std::string(argv[i]) == "-p" || 
                std::string(argv[i]) == "--prefix") {
            i++;
            if (i == argc) {
                ERROR("Option '-p' used without providing prefix for storing results");
            }
            try {
                if (argv[i][0] == '/' || std::string(argv[i]).substr(0, 2) == "./" ||
                        std::string(argv[i]).substr(0, 3) == "../")
                    this->out_prefix = std::string(argv[i++]);
                else
                    this->out_prefix = "./" + std::string(argv[i++]);
                std::string dir = parent_path(out_prefix);
                create_directory(dir);
            } catch (const std::exception & e) {
                ERROR("%s", e.what());
            }
/**************************************************************************************************/
        } else if (std::string(argv[i]) == "-f" || 
                std::string(argv[i]) == "--filter") {
            i++;
            if (i == argc) {
                ERROR("Option '--filter' used without providing filters");
            }
            try {
                size_t filters_before = this->filters.size();
                std::stringstream filters_ss(argv[i++]);
                std::string filter;
                while (getline(filters_ss, filter, ',')) {
                    // leading, interior, and trailing commas leave an empty field, naming nothing
                    if (filter.empty()) continue;
                    this->filters.push_back(filter);
                    this->filter_ids.push_back(-1);
                }
                // an empty filter list means "keep everything", the opposite of what was asked for
                if (this->filters.size() == filters_before)
                    ERROR("Option '--filter' provided no filter names");
            } catch (const std::exception & e) {
                ERROR("%s", e.what());
            }
/**************************************************************************************************/
        } else if (std::string(argv[i]) == "-l" || 
                std::string(argv[i]) == "--largest-variant") {
            i++;
            if (i == argc) {
                ERROR("Option '--largest-variant' used without providing maximum variant size");
            }
            try {
                this->max_size = std::stoi(argv[i++]);
            } catch (const std::exception & e) {
                ERROR("Invalid maximum variant size provided");
            }
/**************************************************************************************************/
        } else if (std::string(argv[i]) == "-sv" || 
                std::string(argv[i]) == "--sv-threshold") {
            i++;
            if (i == argc) {
                ERROR("Option '--sv-threshold' used without providing an SV threshold size");
            }
            try {
                this->sv_threshold = std::stoi(argv[i++]);
            } catch (const std::exception & e) {
                ERROR("Invalid SV threshold size provided");
            }
            if (this->sv_threshold < 2) {
                ERROR("Must provide larger SV threshold size");
            }
/**************************************************************************************************/
        } else if (std::string(argv[i]) == "-q" ||
                std::string(argv[i]) == "--min-qual") {
            i++;
            if (i == argc) {
                ERROR("Option '--min-qual' used without providing minimum variant quality");
            }
            try {
                this->min_qual = std::stoi(argv[i++]);
            } catch (const std::exception & e) {
                ERROR("Invalid minimum variant quality provided");
            }
            if (this->min_qual < 0) {
                ERROR("Must provide non-negative minimum variant quality");
            }
/**************************************************************************************************/
        } else if (std::string(argv[i]) == "-mq" ||
                std::string(argv[i]) == "--max-qual") {
            i++;
            if (i == argc) {
                ERROR("Option '--max-qual' used without providing maximum variant quality");
            }
            try {
                this->max_qual = std::stoi(argv[i++]);
            } catch (const std::exception & e) {
                ERROR("Invalid maximum variant quality provided");
            }
/**************************************************************************************************/
        } else if (std::string(argv[i]) == "-n" || 
                std::string(argv[i]) == "--no-output-files") {
            i++;
            this->write = false;
/**************************************************************************************************/
        } else if (std::string(argv[i]) == "-h" || 
                std::string(argv[i]) == "--help") {
            i++;
            print_help = true;
/**************************************************************************************************/
        } else if (std::string(argv[i]) == "--version") {
            i++;
            this->print_version();
/**************************************************************************************************/
        } else if (std::string(argv[i]) == "-x" || 
                std::string(argv[i]) == "--mismatch-penalty") {
            i++;
            if (i == argc) ERROR("Option '-x' used without providing mismatch penalty");
            try {
                this->sub = std::stoi(argv[i++]);
            } catch (const std::exception & e) {
                ERROR("Invalid mismatch penalty provided");
            }
            if (this->sub < 0) {
                ERROR("Must provide non-negative mismatch penalty");
            }
/**************************************************************************************************/
        } else if (std::string(argv[i]) == "-o" || 
                std::string(argv[i]) == "--gap-open-penalty") {
            i++;
            if (i == argc) ERROR("Option '-o' used without providing gap-opening penalty");
            try {
                this->open = std::stoi(argv[i++]);
            } catch (const std::exception & e) {
                ERROR("Invalid gap-opening penalty provided");
            }
            if (this->open < 0) {
                ERROR("Must provide non-negative gap-opening penalty");
            }
/**************************************************************************************************/
        } else if (std::string(argv[i]) == "-e" || 
                std::string(argv[i]) == "--gap-extend-penalty") {
            i++;
            if (i == argc) ERROR("Option '-e' used without providing gap-extension penalty");
            try {
                this->extend = std::stoi(argv[i++]);
            } catch (const std::exception & e) {
                ERROR("Invalid gap-extension penalty provided");
            }
            if (this->extend < 0) {
                ERROR("Must provide non-negative gap-extension penalty");
            }
/**************************************************************************************************/
        } else if (std::string(argv[i]) == "-i" || 
                std::string(argv[i]) == "--max-iterations") {
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
/**************************************************************************************************/
        } else if (std::string(argv[i]) == "-s" ||
                std::string(argv[i]) == "--max-supercluster-size") {
            i++;
            if (i == argc) ERROR("Option '-s' used without providing max supercluster size");
            try {
                this->max_supercluster_size = std::stoi(argv[i++]);
            } catch (const std::exception & e) {
                ERROR("Invalid max supercluster size provided");
            }
            if (this->max_supercluster_size < 1) {
                ERROR("Max supercluster size must be positive");
            }
/**************************************************************************************************/
        } else if (std::string(argv[i]) == "-t" ||
                std::string(argv[i]) == "--max-threads") {
            i++;
            if (i == argc) {
                ERROR("Option '--max-threads' used without providing max threads");
            }
            try {
                this->max_threads = std::stoi(argv[i++]);
            } catch (const std::exception & e) {
                ERROR("Invalid max threads provided");
            }
            if (this->max_threads < 1) {
                ERROR("Max threads must be positive");
            }
/**************************************************************************************************/
        } else if (std::string(argv[i]) == "-ct" ||
                std::string(argv[i]) == "--credit-threshold") {
            i++;
            if (i == argc) {
                ERROR("Option '--credit-threshold' used without providing value");
            }
            try {
                this->credit_threshold = std::stod(argv[i++]);
            } catch (const std::exception & e) {
                ERROR("Invalid credit threshold provided");
            }
            if (this->credit_threshold <= 0 || this->credit_threshold > 1) {
                ERROR("Provided credit threshold must be on the interval (0,1]");
            }
/**************************************************************************************************/
        } else if (std::string(argv[i]) == "-r" ||
                std::string(argv[i]) == "--max-ram") {
            i++;
            size_t offset = 0;
            if (i == argc) {
                ERROR("Option '--max-ram' used without providing max RAM");
            }
            try {
                this->max_ram = std::stod(argv[i], &offset);
                if (offset != std::string(argv[i]).size()) {
                    WARN("Trailing text '%s' ignored for --max-ram, units are GB",
                            std::string(argv[i]).substr(offset).data());
                }
                i++;
            } catch (const std::exception & e) {
                ERROR("Invalid max RAM provided");
            }
            if (this->max_ram <= 0) {
                ERROR("Max RAM must be positive");
            }
/**************************************************************************************************/
        } else if (std::string(argv[i]) == "-ci" || 
                std::string(argv[i]) == "--citation") {
            i++;
            print_cite = true;
/**************************************************************************************************/
        } else if (std::string(argv[i]) == "-v" ||
                std::string(argv[i]) == "--verbosity") {
            i += 2; // already handled
        } else {
            ERROR("Unexpected option '%s'", argv[i]);
        }
    }

    // final checks, independent of order command-line params are set
    if (this->max_qual < this->min_qual) {
        ERROR("Maximum variant quality must exceed minimum variant quality");
    }

    // warn about variant size exclusions and categories
    if (this->max_size < this->sv_threshold) {
        WARN("No SVs will be evaluated, since --largest-variant %d < --sv-threshold %d", 
                this->max_size, this->sv_threshold);
    }
    if (this->max_supercluster_size < this->max_size + 2) {
        ERROR("Invalid option selected: --max-supercluster-size %d < --largest-variant %d + 2",
                this->max_supercluster_size, this->max_size);
    }
    if (this->max_size == 1) {
        WARN("Only SNPs will be evaluated with --largest-variant %d", this->max_size);
    }

    // recalculate thread/RAM steps, now that --max-threads and --max-ram are known
    this->set_thread_ram_steps();

    if (print_help)
        this->print_usage();
    else if (print_cite)
        this->print_citation();
}

/* --------------------------------------------------------------------------- */

/**
 * @brief Recomputes the thread and RAM scheduling steps from max_threads and max_ram.
 *
 * Each step halves the thread count, so each doubles the RAM available per thread. The
 * constructor calls this so that the steps are already consistent with the default max_threads
 * and max_ram before parse_args() runs; parse_args() calls it again once the options are known.
 * Both vectors are cleared first, so repeated calls replace the steps rather than appending to
 * them.
 */
void Globals::set_thread_ram_steps()
{
    this->thread_steps.clear();
    this->ram_steps.clear();
    this->thread_nsteps = 0;
    int threads = this->max_threads;
    while (threads > 0) {
        this->thread_steps.push_back(threads);
        this->ram_steps.push_back(this->max_ram / threads);
        this->thread_nsteps++;
        threads /= 2;
    }
}

/* --------------------------------------------------------------------------- */

/** @brief Prints program version string to stdout. */
void Globals::print_version() const
{
    printf("%s v%s\n", this->PROGRAM.data(), this->VERSION.data());
}

/** @brief Prints usage information and all command-line options to stdout. */
void Globals::print_usage() const
{
    printf("Usage: vcfdist <query.vcf> <truth.vcf> <ref.fasta> [options]\n"); 

    printf("\nRequired:\n");
    printf("  <STRING>\tquery.vcf\tphased VCF file containing variant calls to evaluate \n");
    printf("  <STRING>\ttruth.vcf\tphased VCF file containing ground truth variant calls \n");
    printf("  <STRING>\tref.fasta\tFASTA file containing draft reference sequence \n");

    printf("\nOptions:\n");
    printf("\n  Inputs/Outputs:\n");
    printf("  -b, --bed <STRING>\n");
    printf("      BED file containing regions to evaluate\n");
    printf("  -v, --verbosity <INTEGER> [%d]\n", this->verbosity);
    printf("      printing verbosity (0: succinct, 1: default, 2:verbose)\n");
    printf("  -p, --prefix <STRING> [./]\n");
    printf("      prefix for output files (directories need a trailing slash)\n");
    printf("  -n, --no-output-files\n");
    printf("      skip writing output files, only print summary to console\n");

    printf("\n  Variant Filtering/Selection:\n");
    printf("  -f, --filter <STRING1,STRING2...> [ALL]\n");
    printf("      select just variants with these FILTER values (OR operation)\n");
    printf("  -l, --largest-variant <INTEGER> [%d]\n", this->max_size);
    printf("      maximum variant size; larger variants are ignored\n");
    printf("  -sv, --sv-threshold <INTEGER> [%d]\n", this->sv_threshold);
    printf("      variants of this size or larger are considered SVs, not INDELs\n");
    printf("  -q, --min-qual <INTEGER> [%d]\n", this->min_qual);
    printf("      minimum variant quality; lower quality variants are ignored\n");
    printf("  -mq, --max-qual <INTEGER> [%d]\n", this->max_qual);
    printf("      maximum variant quality; higher variant qualities are thresholded\n");

    printf("\n  Clustering:\n");
    printf("  -s, --max-supercluster-size <INTEGER> [%d]\n", this->max_supercluster_size);
    printf("      maximum supercluster size (larger superclusters are split)\n");
    /* printf("  -i, --max-iterations <INTEGER> [%d]\n", g.max_cluster_itrs); */
    /* printf("      maximum iterations for expanding/merging clusters\n"); */
    /* printf("  -x, --mismatch-penalty <INTEGER> [%d]\n", g.sub); */
    /* printf("      Smith-Waterman mismatch (substitution) penalty\n"); */
    /* printf("  -o, --gap-open-penalty <INTEGER> [%d]\n", g.open); */
    /* printf("      Smith-Waterman gap opening penalty\n"); */
    /* printf("  -e, --gap-extend-penalty <INTEGER> [%d]\n", g.extend); */
    /* printf("      Smith-Waterman gap extension penalty\n"); */

    printf("\n  Precision-Recall:\n");
    printf("  -ct, --credit-threshold <FLOAT> [%.2f]\n", this->credit_threshold);
    printf("      minimum partial credit to consider a variant a true positive\n");

    printf("\n  Resource Usage:\n");
    printf("  -t, --max-threads <INTEGER> [%d]\n", this->max_threads);
    printf("      maximum threads to use for clustering and precision/recall alignment\n");
    printf("  -r, --max-ram <FLOAT> [%.2fGB]\n", this->max_ram);
    printf("      (approximate) maximum RAM to use for precision/recall alignment\n");

    printf("\n  Miscellaneous:\n");
    printf("  -h, --help\n");
    printf("      show this help message\n");
    printf("  -ci, --citation\n");
    printf("      please cite vcfdist if used in your analyses. Thanks :)\n");
    printf("  -v, --version\n");
    printf("      print %s version (v%s)\n", this->PROGRAM.data(), this->VERSION.data());
}



/**
 * @brief Initializes one named timer object per pipeline stage.
 * @param[in] timer_strs Vector of timer names matching TIME_* constant indices
 */
void Globals::init_timers(const std::vector<std::string> & timer_strs) {
    for (const std::string & timer_name : timer_strs) {
        this->timers.push_back( timer(timer_name) );
    }
}


/**
 * @brief Returns the timer for one pipeline stage.
 * @param[in] t Pipeline stage
 * @return Reference to that stage's timer
 */
timer & Globals::stage(timer_t t) {
    return this->timers[idx(t)];
}


/** @brief Prints publication citation in MLA and BibTeX formats. */
void Globals::print_citation() const
{
    printf("\nMLA Format:\n\n");
    printf("Dunn, Tim, et al. \"Jointly benchmarking small and structural variant calls with vcfdist.\" Genome Biology 25.1 (2024): 253.\n");
    printf("\nBibTeX Format:\n\n");
    printf("@article{dunn2024vcfdist,\n");
    printf("  title={Jointly benchmarking small and structural variant calls with vcfdist},\n");
    printf("  author={Dunn, Tim and Zook, Justin M and Holt, James M and Narayanasamy, Satish},\n");
    printf("  journal={Genome Biology},\n");
    printf("  year={2024},\n");
    printf("  volume={25},\n");
    printf("  number={1},\n");
    printf("  pages={253},\n");
    printf("  publisher={Springer},\n");
    printf("  doi={10.1186/s13059-024-03394-5},\n");
    printf("  URL={https://doi.org/10.1186/s13059-024-03394-5}\n");
    printf("}\n");
}


/**
 * @brief Extracts the parent directory path from a file path string.
 * @param[in] out_prefix Full file path or prefix string
 * @return Parent directory path with trailing slash, or empty string if no parent
 */
std::string parent_path(const std::string & out_prefix) {
    for (int i = out_prefix.size()-1; i >= 0; i--) {
        if (out_prefix[i] == '/')
            return out_prefix.substr(0, i+1);
    }
    return "";
}


/**
 * @brief Creates a directory and all necessary parent directories.
 * @param[in] dir Full directory path to create; an empty path names nothing and is a no-op
 * @throws Errors if directory creation fails for reasons other than EEXIST
 */
void create_directory(const std::string & dir) {
    // an empty path would leave the scan below starting one byte past the end of the copy
    if (dir.empty()) return;

    char *p = strdup(dir.data());
    char *sep = strchr(p+1, '/');
    while(sep != NULL) {
        *sep = '\0';
        if (mkdir(p, 0755) && errno != EEXIST) {
            ERROR("Unable to create directory '%s'", p);
        }
        *sep = '/';
        sep = strchr(sep+1, '/');
    }
    free(p);
}
