#include "edit.h"
#include "globals.h"

void editData::add_edits(const std::string & ctg, int pos, uint8_t hap,
        const std::vector<int> & cig, int sc, int qual) {
   int cig_ptr = 0;
   int type = PTR_DIAG;
   int last_type = PTR_DIAG;
   int len = 0;
   int var_pos = pos;

   while (cig_ptr < int(cig.size())) {
       type = cig[cig_ptr];
       
       if (type != last_type) {
           
           // save last type
           switch (last_type) {
               case PTR_DIAG: // do nothing
                   break;
               case PTR_SUB:
                   add_edit(ctg, var_pos-1, hap, TYPE_SUB, 1, sc, qual);
                   break;
               case PTR_INS:
                   add_edit(ctg, var_pos, hap, TYPE_INS, len, sc, qual);
                   break;
               case PTR_DEL:
                   add_edit(ctg, var_pos-len, hap, TYPE_DEL, len, sc, qual);
                   break;
           }

           // update pointers
           switch (type) {
               case PTR_DIAG:
                   var_pos++;
                   cig_ptr += 2;
                   break;
               case PTR_SUB:
                   var_pos++;
                   cig_ptr += 2;
                   break;
               case PTR_INS:
                   cig_ptr++;
                   break;
               case PTR_DEL:
                   var_pos++;
                   cig_ptr++;
                   break;
           }
           
           last_type = type;
           len = 1;

       } else { // type == last_type
           switch (type) {
               case PTR_DIAG:
                   var_pos++;
                   cig_ptr += 2;
                   break;
               case PTR_SUB:
                   add_edit(ctg, var_pos-1, hap, TYPE_SUB, 1, sc, qual);
                   var_pos++;
                   cig_ptr += 2;
                   break;
               case PTR_INS:
                   cig_ptr++;
                   break;
               case PTR_DEL:
                   var_pos++;
                   cig_ptr++;
                   break;
           }
           len++;
       }

   }
}

void editData::add_edit(const std::string & ctg, int pos, uint8_t hap,
        uint8_t type, int len, int sc, int qual) {
    this->ctgs.push_back(ctg);
    this->poss.push_back(pos);
    this->lens.push_back(len);
    this->haps.push_back(hap);
    this->types.push_back(type);
    this->superclusters.push_back(sc);
    this->quals.push_back(qual);
    this->n++;
}

/*******************************************************************************/

int editData::get_ed(int qual, int type) const {
    int edit_dist = 0;
    for (int i = 0; i < n; i++) {
        if (this->quals[i] == qual && 
                (type == TYPE_ALL || type == this->types[i])) {
            edit_dist += this->lens[i];
        }
    }
    return edit_dist;
}

int editData::get_de(int qual, int type) const {
    int distinct_edits = 0;
    for (int i = 0; i < n; i++) {
        if (this->quals[i] == qual && 
                (type == TYPE_ALL || type == this->types[i])) {
            distinct_edits++;
        }
    }
    return distinct_edits;
}

int editData::get_score(int qual) const {
    int score = 0;
    for (int i = 0; i < n; i++) {
        if (this->quals[i] == qual) {
            switch (this->types[i]) {
                case TYPE_SUB:
                    score += g.eval_sub;
                    break;
                case TYPE_INS:
                case TYPE_DEL:
                    score += g.eval_open + g.eval_extend*this->lens[i];
                    break;
            }
        }
    }
    return score;
}
