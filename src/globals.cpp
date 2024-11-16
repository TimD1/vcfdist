#include "htslib/vcf.h"

#include <sstream>
#include <sys/stat.h>

#include "globals.h"
#include "print.h"
#include "timer.h"

Globals g;
std::vector<std::string> aln_strs = 
    {"QUERY1-TRUTH1", "QUERY1-TRUTH2", "QUERY2-TRUTH1", "QUERY2-TRUTH2"};
std::vector<std::string> callset_strs = {"QUERY", "TRUTH"};
std::vector<std::string> error_strs = {"TP", "FP", "FN", "PE", "GE", "??"};
std::vector<std::string> gt_strs = 
    {"0", "1", "0|0", "0|1", "1|0", "1|1", "1|2", "2|1", ".|.", "M|N" };
std::vector<std::string> phase_strs = {"0", "1", "."};
std::vector<std::string> ac_strs = {".", "+", "-"};
std::vector<std::string> region_strs = {"OUTSIDE", "INSIDE ", "BORDER ", "OFF CTG"};
std::vector<std::string> switch_strs = 
    {"FLIP", "SWITCH", "SWITCH+FLIP", "SWITCH_ERR", "FLIP_BEG", "FLIP_END", "NONE"};
std::vector<std::string> timer_strs = 
    {"reading", "clustering", "realigning", "reclustering", "superclustering", "precision/recall", 
     "phasing", "writing", "total"};
std::vector<std::string> type_strs = {"REF", "SNP", "INS", "DEL", "CPX"};
std::vector<std::string> type_strs2 = {"ALL", "SNP", "INS", "DEL", "INDEL"};
std::vector<std::string> vartype_strs = {"SNP", "INDEL", "SV", "ALL"};
 
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

    g.cmd = std::string(argv[0]);
    for (int i = 1; i < argc; i++)
        g.cmd += " " + std::string(argv[i]);
    if (g.verbosity >= 1) INFO("Command: '%s'", g.cmd.data());

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
    if (g.verbosity >= 1) {
        INFO(" ");
        INFO("%s[0/7] Loading reference FASTA%s '%s'", COLOR_PURPLE,
                COLOR_WHITE, ref_fasta_fn.data());
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
/*******************************************************************************/
        } else if (std::string(argv[i]) == "-f" || 
                std::string(argv[i]) == "--filter") {
            i++;
            if (i == argc) {
                ERROR("Option '--filter' used without providing filters");
            }
            try {
                std::stringstream filters_ss(argv[i++]);
                while (filters_ss.good()) {
                    std::string filter;
                    getline(filters_ss, filter, ',');
                    this->filters.push_back(filter);
                    this->filter_ids.push_back(-1);
                }
            } catch (const std::exception & e) {
                ERROR("%s", e.what());
            }
/*******************************************************************************/
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
/*******************************************************************************/
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
            if (g.sv_threshold < 2) {
                ERROR("Must provide larger SV threshold size");
            }
/*******************************************************************************/
        } else if (std::string(argv[i]) == "-mn" ||
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
            if (g.min_qual < 0) {
                ERROR("Must provide non-negative minimum variant quality");
            }
/*******************************************************************************/
        } else if (std::string(argv[i]) == "-mx" ||
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
/*******************************************************************************/
        } else if (std::string(argv[i]) == "-n" || 
                std::string(argv[i]) == "--no-output-files") {
            i++;
            g.write = false;
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
        } else if (std::string(argv[i]) == "-c" || 
                std::string(argv[i]) == "--cluster") {
            i++;
            if (i == argc) {
                ERROR("Option '-c' used without providing clustering method.");
            }
            if (std::string(argv[i]) == "biwfa") {
                g.cluster_method = argv[i];
                i++;
            } else if (std::string(argv[i]) == "size" || 
                    std::string(argv[i]) == "gap") {
                g.cluster_method = argv[i];
                i++;
                if (i == argc) {
                    ERROR("Option '--cluster %s' used without providing minimum gap between (super)clusters", g.cluster_method.data());
                }
                try {
                    this->cluster_min_gap = std::stoi(argv[i++]);
                } catch (const std::exception & e) {
                    ERROR("Invalid minimum gap between (super)clusters provided");
                }
                if (g.cluster_min_gap <= 0) {
                    ERROR("Must provide positive minimum gap between (super)clusters");
                }
            } else {
                ERROR("Invalid clustering option '%s' provided, must be one of: biwfa, size, gap", argv[i]);
            }

/*******************************************************************************/
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
/*******************************************************************************/
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
/*******************************************************************************/
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
/*******************************************************************************/
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
/*******************************************************************************/
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
/*******************************************************************************/
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
/*******************************************************************************/
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
            if (this->max_ram < 0) {
                ERROR("Max RAM must be positive");
            }
/*******************************************************************************/
        } else if (std::string(argv[i]) == "-ro" || 
                std::string(argv[i]) == "--realign-only") {
            i++;
            g.realign_only = true;
/*******************************************************************************/
        } else if (std::string(argv[i]) == "-rt" ||
                std::string(argv[i]) == "--realign-truth") {
            i++;
            g.realign_truth = true;
/*******************************************************************************/
        } else if (std::string(argv[i]) == "-rq" ||
                std::string(argv[i]) == "--realign-query") {
            i++;
            g.realign_query = true;
/*******************************************************************************/
        } else if (std::string(argv[i]) == "-ci" || 
                std::string(argv[i]) == "--citation") {
            i++;
            print_cite = true;
/*******************************************************************************/
        } else if (std::string(argv[i]) == "-v" ||
                std::string(argv[i]) == "--verbosity") {
            i += 2; // already handled
        } else {
            ERROR("Unexpected option '%s'", argv[i]);
        }
    }

    // final checks, independent of order command-line params are set
    if (g.max_qual < g.min_qual) {
        ERROR("Maximum variant quality must exceed minimum variant quality");
    }

    // warn about variant size exclusions and categories
    if (g.max_size < g.sv_threshold) {
        WARN("No SVs will be evaluated, since --largest %d < --sv-threshold %d", 
                g.max_size, g.sv_threshold);
    }
    if (g.max_size == 1) {
        WARN("Only SNPs will be evaluated with --largest %d", g.max_size);
    }

    // calculate thread/RAM steps
    thread_nsteps = 0;
    int threads = this->max_threads;
    while (threads > 0) {
        thread_steps.push_back(threads);
        ram_steps.push_back(this->max_ram / threads);
        thread_nsteps++;
        threads /= 2;
    }

    if (print_help) 
        this->print_usage();
    else if (print_cite)
        this->print_citation();
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
    printf("  <STRING>\tquery.vcf\tphased VCF file containing variant calls to evaluate \n");
    printf("  <STRING>\ttruth.vcf\tphased VCF file containing ground truth variant calls \n");
    printf("  <STRING>\tref.fasta\tFASTA file containing draft reference sequence \n");

    printf("\nOptions:\n");
    printf("\n  Inputs/Outputs:\n");
    printf("  -b, --bed <STRING>\n");
    printf("      BED file containing regions to evaluate\n");
    printf("  -v, --verbosity <INTEGER> [%d]\n", g.verbosity);
    printf("      printing verbosity (0: succinct, 1: default, 2:verbose)\n");
    printf("  -p, --prefix <STRING> [./]\n");
    printf("      prefix for output files (directory needs a trailing slash)\n");
    printf("  -n, --no-output-files\n");
    printf("      skip writing output files, only print summary to console\n");

    printf("\n  Variant Filtering/Selection:\n");
    printf("  -f, --filter <STRING1,STRING2...> [ALL]\n");
    printf("      select just variants with these FILTER values (OR operation)\n");
    printf("  -l, --largest-variant <INTEGER> [%d]\n", g.max_size);
    printf("      maximum variant size, larger variants ignored\n");
    printf("  -sv, --sv-threshold <INTEGER> [%d]\n", g.sv_threshold);
    printf("      variants of this size or larger are considered SVs, not INDELs\n");
    printf("  -mn, --min-qual <INTEGER> [%d]\n", g.min_qual);
    printf("      minimum variant quality, lower qualities ignored\n");
    printf("  -mx, --max-qual <INTEGER> [%d]\n", g.max_qual);
    printf("      maximum variant quality, higher qualities kept but thresholded\n");

    printf("\n  Re-Alignment:\n");
    printf("  -rq, --realign-query\n");
    printf("      realign query variants using Smith-Waterman parameters\n");
    printf("  -rt, --realign-truth\n");
    printf("      realign truth variants using Smith-Waterman parameters\n");
    printf("  -ro, --realign-only\n");
    printf("      standardize truth and query variant representations, then exit\n");
    printf("  -x, --mismatch-penalty <INTEGER> [%d]\n", g.sub);
    printf("      Smith-Waterman mismatch (substitution) penalty\n");
    printf("  -o, --gap-open-penalty <INTEGER> [%d]\n", g.open);
    printf("      Smith-Waterman gap opening penalty\n");
    printf("  -e, --gap-extend-penalty <INTEGER> [%d]\n", g.extend);
    printf("      Smith-Waterman gap extension penalty\n");

    printf("\n  Clustering:\n");
    printf("  -i, --max-iterations <INTEGER> [%d]\n", g.max_cluster_itrs);
    printf("      maximum iterations for expanding/merging clusters\n");
    printf("  -c, --cluster (biwfa | size <INTEGER> | gap <INTEGER>) [biwfa]\n");
    printf("      select clustering method (see Github Wiki for details)\n");

    printf("\n  Precision-Recall:\n");
    printf("  -ct, --credit-threshold <FLOAT> [%.2f]\n", g.credit_threshold);
    printf("      minimum partial credit to consider variant a true positive\n");

    printf("\n  Utilization:\n");
    printf("  -t, --max-threads <INTEGER> [%d]\n", g.max_threads);
    printf("      maximum threads to use for clustering and precision/recall alignment\n");
    printf("  -r, --max-ram <FLOAT> [%.3fGB]\n", g.max_ram);
    printf("      (approximate) maximum RAM to use for precision/recall alignment\n");

    printf("\n  Miscellaneous:\n");
    printf("  -h, --help\n");
    printf("      show this help message\n");
    printf("  -ci, --citation\n");
    printf("      please cite vcfdist if used in your analyses. Thanks :)\n");
    printf("  -v, --version\n");
    printf("      print %s version (v%s)\n", this->PROGRAM.data(), this->VERSION.data());
}



void Globals::init_timers(const std::vector<std::string> & timer_strs) {
    for (const std::string & timer_name : timer_strs) {
        g.timers.push_back( timer(timer_name) );
    }
}


void Globals::print_citation() const
{
    printf("\nMLA Format:\n\n");
    printf("Dunn, T., Narayanasamy, S. \"vcfdist: accurately benchmarking phased small variant calls in human genomes.\" Nature Communications 14, 8149 (2023). https://doi.org/10.1038/s41467-023-43876-x\n");
    printf("\nBibTeX Format:\n\n");
    printf("@article{dunn2023vcfdist,\n");
    printf("  author={Dunn, Tim and Narayanasamy, Satish},\n");
    printf("  title={vcfdist: Accurately benchmarking phased small variant calls in human genomes},\n");
    printf("  journal={Nature Communications},\n");
    printf("  year={2023},\n");
    printf("  volume={14},\n");
    printf("  number={1},\n");
    printf("  pages={8149},\n");
    printf("  issn={2041-1723},\n");
    printf("  doi={10.1038/s41467-023-43876-x},\n");
    printf("  URL={https://doi.org/10.1038/s41467-023-43876-x}\n");
    printf("}\n");
}


std::string parent_path(const std::string & out_prefix) {
    for (int i = out_prefix.size()-1; i >= 0; i--) {
        if (out_prefix[i] == '/')
            return out_prefix.substr(0, i+1);
    }
    return "";
}


void create_directory(const std::string & dir) {
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
