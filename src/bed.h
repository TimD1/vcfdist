#ifndef _BED_H_
#define _BED_H_

#include <string>
#include <vector>
#include <unordered_map>
#include <fstream>
#include <sstream>
#include <algorithm>

#define OUTSIDE 0
#define INSIDE  1
#define BORDER  2
#define OFF_CTG 3

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

private:
    std::unordered_map<std::string, contigRegions> regions;
    std::vector<std::string> contigs;
};

#endif
