#include "htslib/vcf.h"

#include "globals.h"
#include "print.h"
#include "timer.h"

void Globals::parse_args(int argc, char ** argv) {

    /* if required arguments are not provided, you can only print help and exit */
    bool print_cite = false;
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
            } else if (std::string(argv[i]) == "-v" ||
                    std::string(argv[i]) == "--version") {
                i++;
                this->print_version();
            } else if (std::string(argv[i]) == "-c" || 
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

    for (int i = 0; i < argc; i++)
        g.cmd += " " + std::string(argv[i]);
    if (g.verbosity >= 1) INFO("Command: '%s'", g.cmd.data()+1);

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
        INFO("%s[0/8] Loading reference FASTA%s '%s'", COLOR_PURPLE,
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
                this->out_prefix = std::string(argv[i++]);
                std::filesystem::create_directory(
                        std::filesystem::path(out_prefix).parent_path());
            } catch (const std::exception & e) {
                ERROR("%s", e.what());
            }
/*******************************************************************************/
        } else if (std::string(argv[i]) == "-s" || 
                std::string(argv[i]) == "--smallest-variant") {
            i++;
            if (i == argc) {
                ERROR("Option '--smallest-variant' used without providing minimum variant size");
            }
            try {
                this->min_size = std::stoi(argv[i++]);
            } catch (const std::exception & e) {
                ERROR("Invalid minimum variant size provided");
            }
            if (g.min_size < 0) {
                ERROR("Must provide non-negative minimum variant size");
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
        } else if (std::string(argv[i]) == "--min-qual") {
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
        } else if (std::string(argv[i]) == "--max-qual") {
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
        } else if (std::string(argv[i]) == "-qx" || 
                std::string(argv[i]) == "--query-mismatch-penalty") {
            i++;
            if (i == argc) {
                ERROR("Option '-qx' used without providing mismatch penalty");
            } else if (this->query_penalties_set[PEN_SUB]) {
                ERROR("Query mismatch penalty already set, cannot use '-qx'");
            }
            try {
                this->query_sub = std::stoi(argv[i++]);
                this->query_penalties_set[PEN_SUB] = true;
            } catch (const std::exception & e) {
                ERROR("Invalid query mismatch penalty provided");
            }
            if (this->query_sub < 0) {
                ERROR("Must provide non-negative query mismatch penalty");
            }
/*******************************************************************************/
        } else if (std::string(argv[i]) == "-tx" || 
                std::string(argv[i]) == "--truth-mismatch-penalty") {
            i++;
            if (i == argc) {
                ERROR("Option '-tx' used without providing mismatch penalty");
            } else if (this->truth_penalties_set[PEN_SUB]) {
                ERROR("Truth mismatch penalty already set, cannot use '-tx'");
            }
            try {
                this->truth_sub = std::stoi(argv[i++]);
                this->truth_penalties_set[PEN_SUB] = true;
            } catch (const std::exception & e) {
                ERROR("Invalid truth mismatch penalty provided");
            }
            if (this->truth_sub < 0) {
                ERROR("Must provide non-negative truth mismatch penalty");
            }
/*******************************************************************************/
        } else if (std::string(argv[i]) == "-ex" || 
                std::string(argv[i]) == "--eval-mismatch-penalty") {
            i++;
            if (i == argc) {
                ERROR("Option '-ex' used without providing mismatch penalty");
            } else if (this->eval_penalties_set[PEN_SUB]) {
                ERROR("Evaluation mismatch penalty already set, cannot use '-ex'");
            }
            try {
                this->eval_sub = std::stoi(argv[i++]);
                this->eval_penalties_set[PEN_SUB] = true;
            } catch (const std::exception & e) {
                ERROR("Invalid evaluation mismatch penalty provided");
            }
            if (this->eval_sub < 0) {
                ERROR("Must provide non-negative evaluation mismatch penalty");
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
        } else if (std::string(argv[i]) == "--max-threads") {
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
        } else if (std::string(argv[i]) == "--max-ram") {
            i++;
            if (i == argc) {
                ERROR("Option '--max-ram' used without providing max RAM");
            }
            try {
                this->max_ram = std::stod(argv[i++]);
            } catch (const std::exception & e) {
                ERROR("Invalid max RAM provided");
            }
            if (this->max_ram < 0) {
                ERROR("Max RAM must be positive");
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
        } else if (std::string(argv[i]) == "-t" ||
                std::string(argv[i]) == "--keep-truth") {
            i++;
            g.keep_truth = true;
/*******************************************************************************/
        } else if (std::string(argv[i]) == "-q" ||
                std::string(argv[i]) == "--keep-query") {
            i++;
            g.keep_query = true;
/*******************************************************************************/
        } else if (std::string(argv[i]) == "--simple-cluster") {
            i++;
            g.simple_cluster = true;
/*******************************************************************************/
        } else if (std::string(argv[i]) == "-c" || 
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

    // final check, independent of order command-line params are set
    if (g.max_size < g.min_size) {
        ERROR("Maximum variant size must exceed smallest variant size");
    }
    if (g.max_qual < g.min_qual) {
        ERROR("Maximum variant quality must exceed minimum variant quality");
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
    printf("  -b, --bed <STRING>\n");
    printf("      BED file containing regions to evaluate\n\n");

    printf("  -p, --prefix <STRING> [./]\n");
    printf("      prefix for output files (directory needs a trailing slash)\n\n");

    printf("  -v, --verbosity <INTEGER> [%d]\n", g.verbosity);
    printf("      printing verbosity (0: succinct, 1: default, 2:verbose)\n\n");

    printf("  -r, --realign-only\n");
    printf("      standardize truth and query variant representations, then exit\n\n");

    printf("  -q, --keep-query\n");
    printf("      do not realign query variants, keep original representation\n\n");

    printf("  -t, --keep-truth\n");
    printf("      do not realign truth variants, keep original representation\n\n");

    printf("  -x, --mismatch-penalty <INTEGER> [%d]\n", g.eval_sub);
    printf("      Smith-Waterman mismatch (substitution) penalty\n\n");

    printf("  -o, --gap-open-penalty <INTEGER> [%d]\n", g.eval_open);
    printf("      Smith-Waterman gap opening penalty\n\n");

    printf("  -e, --gap-extend-penalty <INTEGER> [%d]\n", g.eval_extend);
    printf("      Smith-Waterman gap extension penalty\n\n");

    printf("  --min-qual <INTEGER> [%d]\n", g.min_qual);
    printf("      minimum variant quality, lower qualities ignored\n\n");

    printf("  --max-qual <INTEGER> [%d]\n", g.max_qual);
    printf("      maximum variant quality, higher qualities kept but thresholded\n\n");

    printf("  -s, --smallest-variant <INTEGER> [%d]\n", g.min_size);
    printf("      minimum variant size, smaller variants ignored (SNPs are size 1)\n\n");

    printf("  -l, --largest-variant <INTEGER> [%d]\n", g.max_size);
    printf("      maximum variant size, larger variants ignored\n\n");

    printf("  -i, --max-iterations <INTEGER> [%d]\n", g.max_cluster_itrs);
    printf("      maximum iterations for expanding/merging clusters\n\n");

    printf("  -g, --supercluster-gap <INTEGER> [%d]\n", g.cluster_min_gap);
    printf("      minimum base gap between independent superclusters\n\n");

    printf("  --max-threads <INTEGER> [%d]\n", g.max_threads);
    printf("      maximum threads to use for precision/recall alignment\n");
    printf("      (haps*contigs used for wavefront clustering)\n\n");

    printf("  --max-ram <FLOAT> [%.3fGB]\n", g.max_ram);
    printf("      maximum RAM to use for precision/recall alignment\n");
    printf("      (work in-progress, more may be used in other steps)\n\n");

    printf("  -h, --help\n");
    printf("      show this help message\n\n");

    printf("  -a, --advanced\n");
    printf("      show advanced options\n\n");

    printf("  -c, --citation\n");
    printf("      please cite vcfdist if used in your analyses\n\n");

    printf("  -v, --version\n");
    printf("      print %s version (v%s)\n\n", 
            this->PROGRAM.data(), this->VERSION.data());

    if (!this->advanced) return;
    printf("\nAdvanced Options: (not recommended, only for evaluation)\n");
    printf("  --simple-cluster\n");
    printf("      instead of Smith-Waterman clustering, use gap-based clustering \n\n");

    printf("  -qx, --query-mismatch-penalty <INTEGER> [%d]\n", g.query_sub);
    printf("      mismatch penalty (query variant realignment)\n\n");

    printf("  -qo, --query-gap-open-penalty <INTEGER> [%d]\n", g.query_open);
    printf("      gap opening penalty (query variant realignment)\n\n");

    printf("  -qe, --query-gap-extend-penalty <INTEGER> [%d]\n", g.query_extend);
    printf("      gap extension penalty (query variant realignment)\n\n");

    printf("  -tx, --truth-mismatch-penalty <INTEGER> [%d]\n", g.truth_sub);
    printf("      mismatch penalty (truth variant realignment)\n\n");

    printf("  -to, --truth-gap-open-penalty <INTEGER> [%d]\n", g.truth_open);
    printf("      gap opening penalty (truth variant realignment)\n\n");

    printf("  -te, --truth-gap-extend-penalty <INTEGER> [%d]\n", g.truth_extend);
    printf("      gap extension penalty (truth variant realignment)\n\n");

    printf("  -ex, --eval-mismatch-penalty <INTEGER> [%d]\n", g.eval_sub);
    printf("      mismatch penalty (distance evaluation)\n\n");

    printf("  -eo, --eval-gap-open-penalty <INTEGER> [%d]\n", g.eval_open);
    printf("      gap opening penalty (distance evaluation)\n\n");

    printf("  -ee, --eval-gap-extend-penalty <INTEGER> [%d]\n", g.eval_extend);
    printf("      gap extension penalty (distance evaluation)\n\n");

}



void Globals::init_timers(std::vector<std::string> timer_strs) {
    for (std::string timer_name : timer_strs) {
        g.timers.push_back( timer(timer_name) );
    }
}


void Globals::print_citation() const
{
    printf("\nMLA Format:\n\n");
    printf("  Dunn, Tim, and Satish Narayanasamy. \"vcfdist: Accurately benchmarking phased small variant calls in human genomes.\" bioRxiv (2023): 2023-03.\n");
    printf("\nBibTeX Format:\n\n");
    printf("  @article {dunn2023vcfdist,\n");
    printf("    author = {Dunn, Tim and Narayanasamy, Satish},\n");
    printf("    title = {vcfdist: Accurately benchmarking phased small variant calls in human genomes},\n");
    printf("    elocation-id = {2023.03.10.532078},\n");
    printf("    year = {2023},\n");
    printf("    doi = {10.1101/2023.03.10.532078},\n");
    printf("    publisher = {Cold Spring Harbor Laboratory},\n");
    printf("    URL = {https://www.biorxiv.org/content/early/2023/03/12/2023.03.10.532078},\n");
    printf("    eprint = {https://biorxiv.org/content/early/2023/03/12/2023.03.10.532078.full.pdf},\n");
    printf("    journal = {bioRxiv}\n");
    printf("  }\n");
}
