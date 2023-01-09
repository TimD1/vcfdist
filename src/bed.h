#ifndef _BED_H_
#define _BED_H_

#include <string>
#include <vector>
#include <unordered_map>
#include <fstream>
#include <sstream>
#include <algorithm>

#include "defs.h"
#include "fasta.h"
#include "variant.h"

struct contigRegions {
    std::vector<int> starts;
    std::vector<int> stops;
    int n;
};

class bedData {
public:

    // constructors
    bedData() {;}
    bedData(const std::string & bed_fn);

    void add(const std::string & contig, const int & start, const int & stop);
    void check();
    int contains(std::string contig, const int & start, const int & stop);

    operator std::string() const;

    long size;

    std::unordered_map<std::string, contigRegions> regions;
    std::vector<std::string> contigs;
};

void check_contigs(
        std::unique_ptr<variantData> & query_ptr,
        std::unique_ptr<variantData> & truth_ptr,
        std::shared_ptr<fastaData> ref_ptr);

#endif
