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


/*******************************************************************************/


void editData::write_distance() {

    // log all distance results
    std::string dist_fn = g.out_prefix + "distance.tsv";
    FILE* out_dists = 0;
    if (g.write) {
        out_dists = fopen(dist_fn.data(), "w");
        if (g.verbosity >= 1) INFO(" ");
        if (g.verbosity >= 1) INFO("  Writing distance results to '%s'", dist_fn.data());
        fprintf(out_dists, "MIN_QUAL\tSUB_DE\tINS_DE\tDEL_DE\tSUB_ED\tINS_ED\tDEL_ED\t"
                "DISTINCT_EDITS\tEDIT_DIST\tALN_SCORE\tALN_QSCORE\n");
    }

    // get original scores / distance (above g.max_qual, no vars applied)
    std::vector<int> orig_edit_dists(TYPES, 0);
    std::vector<int> orig_distinct_edits(TYPES, 0);
    int orig_score = this->get_score(g.max_qual+1);
    for (int type = 0; type < TYPES; type++) {
        orig_edit_dists[type] = this->get_ed(g.max_qual+1, type);
        orig_distinct_edits[type] = this->get_de(g.max_qual+1, type);
    }

    std::vector<double> best_score(TYPES, std::numeric_limits<double>::max());
    std::vector<int> best_qual(TYPES, 0);
    for (int q = g.min_qual; q <= g.max_qual+1; q++) { // all qualities

        // get ED/DE for each Q threshold, for each type
        std::vector<int> edit_dists(TYPES, 0);
        std::vector<int> distinct_edits(TYPES, 0);
        for (int type = 0; type < TYPES; type++) {
            edit_dists[type] = this->get_ed(q, type);
            distinct_edits[type] = this->get_de(q, type);

            // save best Q threshold so far
            double score = double(edit_dists[type]) * distinct_edits[type];
            if (score < best_score[type]) {
                best_score[type] = score;
                best_qual[type] = q;
            }
        }

        if (g.write) fprintf(out_dists, 
                "%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%f\n", 
                q, distinct_edits[TYPE_SUB],
                distinct_edits[TYPE_INS], distinct_edits[TYPE_DEL],
                edit_dists[TYPE_SUB], edit_dists[TYPE_INS], edit_dists[TYPE_DEL],
                distinct_edits[TYPE_ALL], edit_dists[TYPE_ALL],
                this->get_score(q), qscore(double(this->get_score(q))/orig_score));
    }
    if (g.write) fclose(out_dists);

    // summarize distance results
    std::string dist_summ_fn = g.out_prefix + "distance-summary.tsv";
    FILE* dists_summ = 0;
    if (g.write) {
        dists_summ = fopen(dist_summ_fn.data(), "w");
        if (g.verbosity >= 1) 
            INFO("  Writing distance summary to '%s'", dist_summ_fn.data());
        fprintf(dists_summ, "VAR_TYPE\tTHRESHOLD\tMIN_QUAL\tEDIT_DIST\tDISTINCT_EDITS\tED_QSCORE\tDE_QSCORE\tALN_QSCORE\n");
    }
    INFO(" ");
    INFO("%sALIGNMENT DISTANCE SUMMARY%s", COLOR_BLUE, COLOR_WHITE);
    for (int type = 0; type < TYPES; type++) {

        // skip INS/DEL individually unless higher print verbosity
        if (g.verbosity == 0 && type != TYPE_ALL)
            continue;
        if (g.verbosity == 1 && (type == TYPE_INS || type == TYPE_DEL))
            continue;

        INFO(" ");
        if (type == TYPE_ALL) {
            INFO("%sTYPE\tTHRESHOLD\tEDIT_DIST\tDISTINCT_EDITS\tED_QSCORE\tDE_QSCORE\tALN_QSCORE%s",
                    COLOR_BLUE, COLOR_WHITE);
        } else {
            INFO("%sTYPE\tTHRESHOLD\tEDIT_DIST\tDISTINCT_EDITS\tED_QSCORE\tDE_QSCORE%s",
                    COLOR_BLUE, COLOR_WHITE);
        }
        std::vector<int> quals = {g.min_qual, best_qual[type], g.max_qual+1};
        std::vector<std::string> thresholds = {"NONE", "BEST", "REF "};
        for (int i = 0; i < int(quals.size()); i++) {

            // fill out ED/DE for selected quals
            int q = quals[i];
            std::vector<int> edit_dists(TYPES, 0);
            std::vector<int> distinct_edits(TYPES, 0);
            for (int type = 0; type < TYPES; type++) {
                edit_dists[type] =this->get_ed(q, type);
                distinct_edits[type] = this->get_de(q, type);
            }

            float ed_qscore = qscore(double(edit_dists[type]) / orig_edit_dists[type]);
            float de_qscore = qscore(double(distinct_edits[type]) / orig_distinct_edits[type]);
            float all_qscore = type == TYPE_ALL ? qscore(double(this->get_score(q)) / orig_score) : 0;

            // print summary
            if (g.write) fprintf(dists_summ, "%s\t%s\t%d\t%d\t%d\t%f\t%f\t%f\n", 
                    type_strs2[type].data(), thresholds[i].data(),
                    q, edit_dists[type], distinct_edits[type], ed_qscore, de_qscore, all_qscore);
            if (type == TYPE_ALL) {
                INFO("%s%s\t%s Q >= %d\t%-16d%-16d%f\t%f\t%f%s", 
                    COLOR_BLUE, type_strs2[type].data(), thresholds[i].data(),
                    q, edit_dists[type], distinct_edits[type], ed_qscore, de_qscore, all_qscore, COLOR_WHITE);
            } else {
                INFO("%s%s\t%s Q >= %d\t%-16d%-16d%f\t%f%s", 
                    COLOR_BLUE, type_strs2[type].data(), thresholds[i].data(),
                    q, edit_dists[type], distinct_edits[type], ed_qscore, de_qscore, COLOR_WHITE);
            }
        }
    }
    INFO(" ");
    if (g.write) fclose(dists_summ);
}


/*******************************************************************************/


void editData::write_edits() {
    std::string edit_fn = g.out_prefix + "edits.tsv";
    if (g.verbosity >= 1) INFO("  Writing edit results to '%s'", edit_fn.data());
    FILE* out_edits = fopen(edit_fn.data(), "w");
    fprintf(out_edits, "CONTIG\tSTART\tHAP\tTYPE\tSIZE\tSUPERCLUSTER\tMIN_QUAL\tMAX_QUAL\n");
    for (int i = 0; i < this->n; i++) {
        fprintf(out_edits, "%s\t%d\t%d\t%s\t%d\t%d\t%d\t%d\n", 
                this->ctgs[i].data(), this->poss[i], this->haps[i],
                type_strs[this->types[i]].data(), this->lens[i],
                this->superclusters[i], this->min_quals[i], this->max_quals[i]);
    }
    fclose(out_edits);
}
