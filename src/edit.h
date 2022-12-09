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
            const std::vector<int> & cig, int sc, int qual);
    void add_edit(const std::string & ctg, int pos, uint8_t hap, 
            uint8_t type, int len, int sc, int qual);

    int get_ed(int qual, int type=TYPE_ALL) const; // edit distance
    int get_de(int qual, int type=TYPE_ALL) const; // distinct edits
    int get_score(int qual) const;                 // smith-waterman distance

    // data (all of size n)
    std::vector<std::string> ctgs;  // contigs
    std::vector<int> poss;          // variant start positions (0-based)
    std::vector<uint8_t> haps;      // variant haplotypes
    std::vector<uint8_t> types;     // variant type: NONE, SUB, INS, DEL, GRP
    std::vector<int> lens;          // variant lengths
    std::vector<int> superclusters; // variant superclusters
    std::vector<int> quals;         // variant qualities
    int n = 0;
};

#endif
