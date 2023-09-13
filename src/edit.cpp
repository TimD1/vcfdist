#include "edit.h"
#include "globals.h"

void editData::add_edits(const std::string & ctg, int pos, uint8_t hap,
        const std::vector<int> & cig, int sc, int minq, int maxq) {
   int cig_ptr = 0;
   int type = PTR_MAT;
   int last_type = PTR_MAT;
   int len = 0;
   int var_pos = pos;

   while (cig_ptr < int(cig.size())) {
       type = cig[cig_ptr];
       
       if (type != last_type) {
           
           // save last type
           switch (last_type) {
               case PTR_MAT: // do nothing
                   break;
               case PTR_SUB:
                   add_edit(ctg, var_pos-1, hap, TYPE_SUB, 1, sc, minq, maxq);
                   break;
               case PTR_INS:
                   add_edit(ctg, var_pos, hap, TYPE_INS, len, sc, minq, maxq);
                   break;
               case PTR_DEL:
                   add_edit(ctg, var_pos-len, hap, TYPE_DEL, len, sc, minq, maxq);
                   break;
           }

           // update pointers
           switch (type) {
               case PTR_MAT:
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
               case PTR_MAT:
                   var_pos++;
                   cig_ptr += 2;
                   break;
               case PTR_SUB:
                   add_edit(ctg, var_pos-1, hap, TYPE_SUB, 1, sc, minq, maxq);
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
        uint8_t type, int len, int sc, int min_qual, int max_qual) {
    this->ctgs.push_back(ctg);
    this->poss.push_back(pos);
    this->lens.push_back(len);
    this->haps.push_back(hap);
    this->types.push_back(type);
    this->superclusters.push_back(sc);
    this->min_quals.push_back(min_qual);
    this->max_quals.push_back(max_qual);
    this->n++;
}

/*******************************************************************************/

bool is_type(int type, int category) {
    if (type == TYPE_INDEL || type == TYPE_ALL) ERROR("Invalid type in is_type()");
    if (category == TYPE_ALL) return true;
    if (type == category) return true;
    if (category == TYPE_INDEL && (type == TYPE_INS || type == TYPE_DEL)) return true;
    return false;
}

float qscore(double p_error) {
    return std::min(100.0, std::max(0.0, -10 * std::log10(p_error)));
}

int editData::get_ed(int qual, int type) const {
    int edit_dist = 0;
    for (int i = 0; i < this->n; i++) {
        if (qual >= this->min_quals[i] && qual < this->max_quals[i] && 
                is_type(this->types[i], type)) {
            edit_dist += this->lens[i];
        }
    }
    return edit_dist;
}

int editData::get_de(int qual, int type) const {
    int distinct_edits = 0;
    for (int i = 0; i < this->n; i++) {
        if (qual >= this->min_quals[i] && qual < this->max_quals[i] &&
                is_type(this->types[i], type)) {
            distinct_edits++;
        }
    }
    return distinct_edits;
}

int editData::get_score(int qual) const {
    int score = 0;
    for (int i = 0; i < n; i++) {
        if (qual >= this->min_quals[i] && qual < this->max_quals[i]) {
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
