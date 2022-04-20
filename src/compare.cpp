#include <cstdio>
#include <string>
#include <exception>
#include <filesystem>
#include <unordered_map>
#include <vector>
#include <fstream>
#include <iostream>
#include <sstream>

#include "vcf.h"
#include "vcfutils.h"

std::string VERSION = "0.0.1";
void print_version() {
    printf("vcfdist compare v%s\n", VERSION.data());
}

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
    printf("  -h\t\tshow this help message\n");
    printf("  -v\t\tshow version number\n");
}

/* --------------------------------------------------------------------------- */

struct contigRegions {
    std::vector<int> starts;
    std::vector<int> stops;
};

class bedRegions {
public:

    // constructors
    bedRegions() = default;
    bedRegions(std::string bed_fn) {
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

int parse_args(
        int argc, 
        char ** argv, 
        htsFile *& ref_fasta, 
        htsFile *& calls_vcf, 
        htsFile *& truth_vcf, 
        bedRegions *& bed,
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
        printf("[E::compare] failed to open calls_vcf file '%s'.\n", 
                calls_vcf_fn.data());
        return 1;
    }
    truth_vcf_fn = std::string(argv[2]);
    truth_vcf = bcf_open(truth_vcf_fn.data(), "r");
    if (truth_vcf == NULL) {
        printf("[E::compare] failed to open truth_vcf file '%s'.\n", 
                truth_vcf_fn.data());
        return 1;
    }
    ref_fasta_fn = std::string(argv[3]);
    ref_fasta = bcf_open(ref_fasta_fn.data(), "r");
    if (ref_fasta == NULL) {
        printf("[E::compare] failed to open ref_fasta file '%s'.\n", 
                ref_fasta_fn.data());
        return 1;
    }

    /* handle optional arguments */
    for (int i = 4; i < argc;) {
        if (std::string(argv[i]) == "-b") {
            i++;
            if (i == argc) {
                printf("[E::compare] option '-b' used without providing BED.\n");
                return 1;
            }
            try {
                bed = new bedRegions(std::string(argv[i++]));
            } catch (const std::exception & e) {
                printf("[E::compare] %s\n", e.what());
                return 1;
            }
        }
        else if (std::string(argv[i]) == "-o") {
            i++;
            if (i == argc) {
                printf("[E::compare] option '-o' used without providing DIR.\n");
                return 1;
            }
            try {
                out_dir = std::string(argv[i++]);
                std::filesystem::create_directory(out_dir);
            } catch (const std::exception & e) {
                printf("[E::compare] %s\n", e.what());
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
        else {
            printf("[E::compare] unexpected option '%s'.\n", argv[i]);
            return 1;
        }
    }

    return 0;
}

/* --------------------------------------------------------------------------- */
 
int main(int argc, char **argv) {

    // parse and store command-line args
    htsFile *ref_fasta, *calls_vcf, *truth_vcf;
    bedRegions * bed;
    std::string out_dir;
    int exit = parse_args(argc, argv, 
            ref_fasta, calls_vcf, truth_vcf, bed, out_dir);
    if (exit) std::exit(EXIT_SUCCESS);

    std::cout << std::string(*bed) << std::endl;

        
    /* // read header */
    /* bcf_hdr_t *hdr = bcf_hdr_read(inf); */

    // if hdr == NULL
 
    /* // struc for storing each record */
    /* bcf1_t *rec = bcf_init(); */
    /* if (rec == NULL) { */
    /*     return EXIT_FAILURE; */
    /* } */
 
    /* while (bcf_read(inf, hdr, rec) == 0) { */
    /*     if (bcf_is_snp(rec)) { */
    /*         printf("Number of alleles: %lu\n", (unsigned long)rec->n_allele); */
    /*         printf("Chr name is %s\n", bcf_hdr_id2name(hdr, rec->rid)); */
    /*         printf("Position for this SNP is %li\n", rec->pos); */
    /*         printf("REF:%s ALT:%s\n", rec->d.allele[0], rec->d.allele[1]); */
    /*         // iterating over alleles */
    /*         for (int i=0; i<rec->n_allele; ++i) */
    /*             printf("%s\n", rec->d.allele[i]); */
    /*     } else { */
    /*         continue; */
    /*     } */
    /* } */
 
    /* bcf_hdr_destroy(hdr); */
    /* bcf_close(inf); */
    /* bcf_destroy(rec); */
    /* delete bed */

    return EXIT_SUCCESS;
}
