#ifndef _EDIT_H_
#define _EDIT_H_

#include <string>
#include <vector>
#include <cmath>

#include "defs.h"

float qscore(double p_error);

class editData {
public:

    // constructor
    editData() {};

    // helper functions
    void add_edits(const std::string & ctg, int pos, uint8_t hap, 
            const std::vector<int> & cig, int sc, int min_qual, int max_qual);
    void add_edit(const std::string & ctg, int pos, uint8_t hap, 
            uint8_t type, int len, int sc, int min_qual, int max_qual);

    int get_ed(int qual, int type=TYPE_ALL) const; // edit distance
    int get_de(int qual, int type=TYPE_ALL) const; // distinct edits
    int get_score(int qual) const;                 // smith-waterman distance

    // data (all of size n)
    int n = 0;
    std::vector<std::string> ctgs;  // contigs
    std::vector<int> poss;          // edit start position (0-based)
    std::vector<uint8_t> haps;      // edit haplotype
    std::vector<uint8_t> types;     // edit type: SUB, INS, DEL
    std::vector<int> lens;          // edit length
    std::vector<int> superclusters; // edit supercluster
    std::vector<int> min_quals;     // edit quality start range (inclusive)
    std::vector<int> max_quals;     // edit quality end range (exclusive)
};

#endif
