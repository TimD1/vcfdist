#include <filesystem>

#include "htslib/vcf.h"

#include "globals.h"
#include "print.h"

void Globals::parse_args(int argc, char ** argv) {

    /* print help and exit */
    if (argc < 4) {
        if (argc == 2 && std::string(argv[1]) == "-v") {
            this->print_version();
            std::exit(0);
        }
        this->print_usage();
        std::exit(0);
    }

    if (argv[1][0] == '-' || argv[2][0] == '-' || argv[3][0] == '-') {
        WARN("optional arguments should be provided AFTER mandatory arguments");
        this->print_usage();
    }

    /* verify input VCF/FASTA filepaths */
    this->calls_vcf_fn = std::string(argv[1]);
    this->calls_vcf_path = std::filesystem::path(this->calls_vcf_fn);
    htsFile* calls_vcf_fp = bcf_open(calls_vcf_fn.data(), "r");
    if (calls_vcf_fp == NULL) {
        ERROR("failed to open calls_vcf file '%s'", calls_vcf_fn.data());
    } else {
        bcf_close(calls_vcf_fp);
    }

    this->truth_vcf_fn = std::string(argv[2]);
    this->truth_vcf_path = std::filesystem::path(this->truth_vcf_fn);
    htsFile* truth_vcf_fp = bcf_open(truth_vcf_fn.data(), "r");
    if (truth_vcf_fp == NULL) {
        ERROR("failed to open truth_vcf file '%s'", truth_vcf_fn.data());
    } else {
        bcf_close(truth_vcf_fp);
    }

    // load reference FASTA
    this->ref_fasta_fn = std::string(argv[3]);
    INFO("loading reference FASTA '%s'", ref_fasta_fn.data());
    this->ref_fasta_fp = fopen(ref_fasta_fn.data(), "r");
    if (ref_fasta_fp == NULL) {
        ERROR("failed to open ref_fasta file '%s'", ref_fasta_fn.data());
    }

    /* handle optional arguments */
    for (int i = 4; i < argc;) {
        if (std::string(argv[i]) == "-b") {
            i++;
            if (i == argc) {
                ERROR("option '-b' used without providing BED");
            }
            try {
                this->bed_fn = std::string(argv[i]);
                this->bed = bedData(std::string(argv[i++]));
                this->bed.check();
            } catch (const std::exception & e) {
                ERROR("invalid BED file provided");
            }
        }
        else if (std::string(argv[i]) == "-o") {
            i++;
            if (i == argc) {
                ERROR("option '-o' used without providing PREFIX");
            }
            try {
                this->out_prefix = std::string(argv[i++]);
                std::filesystem::create_directory(
                        std::filesystem::path(out_prefix).parent_path());
            } catch (const std::exception & e) {
                ERROR("%s", e.what());
            }
        }
        else if (std::string(argv[i]) == "-q") {
            i++;
            if (i == argc) {
                ERROR("option '-q' used without providing MIN_QUAL");
            }
            try {
                this->min_qual = std::stoi(argv[i++]);
            } catch (const std::exception & e) {
                ERROR("invalid minimum variant quality provided");
            }
        }
        else if (std::string(argv[i]) == "-p") {
            i++;
            if (i == argc) {
                ERROR("option '-p' used without providing VERBOSITY");
            }
            try {
                this->print_verbosity = std::stoi(argv[i++]);
            } catch (const std::exception & e) {
                ERROR("invalid print verbosity provided");
            }
            if (this->print_verbosity < 0 || this->print_verbosity > 2) {
                ERROR("print verbosity %d not a valid option (0,1,2)", 
                        this->print_verbosity);
            }
        }
        else if (std::string(argv[i]) == "-h") {
            i++;
            this->print_usage();
        }
        else if (std::string(argv[i]) == "-v") {
            i++;
            this->print_version();
        }
        else if (std::string(argv[i]) == "-g") {
            i++;
            if (i == argc) {
                ERROR("option '-g' used without providing GAP size");
            }
            try {
                this->gap = std::stoi(argv[i++]);
            } catch (const std::exception & e) {
                ERROR("invalid gap size provided");
            }
            if (g.gap <= 0) {
                ERROR("must provide positive gap size");
            }
        }
        else {
            ERROR("unexpected option '%s'", argv[i]);
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
    printf("Usage: vcfdist <calls.vcf> <truth.vcf> <ref.fasta> [options]\n"); 
    printf("\nRequired:\n");
    printf("  FILE\tcalls.vcf\tphased VCF file containing variant calls to evaluate \n");
    printf("  FILE\ttruth.vcf\tphased VCF file containing ground truth variant calls \n");
    printf("  FILE\tref.fasta\tFASTA file containing reference sequence \n");
    printf("\nOptions:\n");
    printf("  -b FILE\tBED file containing regions to evaluate\n");
    printf("  -o PREFIX\toutput filepath prefix (dir should contain trailing slash)\n");
    printf("  -q MIN_QUAL\tminimum variant quality [0]\n");
    printf("  -g GAP\tminimum base gap between independent clusters [50]\n");
    printf("  -p VERBOSITY\tprinting verbosity (0,1,2) [0]\n");
    printf("  -h\t\tshow this help message\n");
    printf("  -v\t\tshow version number\n");
}
