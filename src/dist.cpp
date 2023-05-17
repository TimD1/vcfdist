#include <string>
#include <set>
#include <vector>
#include <cstdio>
#include <chrono>
#include <utility>
#include <queue>

#include "dist.h"
#include "edit.h"
#include "print.h"
#include "cluster.h"

template <typename T>
bool contains(const std::unordered_set<T> & wave, const T & idx) {
    return wave.find(idx) != wave.end();
}
template <typename T, typename U>
bool contains(const std::unordered_map<T,U> & wave, const T & idx) {
    return wave.find(idx) != wave.end();
}

template <typename S, typename T>
void set(std::unordered_map<S,T> & ptrs, const S & idx, const T & flag) {
    if (ptrs.find(idx) != ptrs.end()) // update in map
        ptrs[idx] |= flag;
    else // add to map
        ptrs[idx] = flag;
}

template <typename S, typename T>
bool test(const std::unordered_map<S,T> & ptrs, const S & idx, const T & flag) {
    return ptrs.find(idx) != ptrs.end() && ptrs.at(idx) & flag;
}


/******************************************************************************/


/* Reverse two strings and their associated pointers for reverse alignment. */
void reverse_ptrs_strs(
        std::string & query, std::string & ref,
        std::vector< std::vector<int> > & query_ptrs, 
        std::vector< std::vector<int> > & ref_ptrs) {

    // reverse strings and vectors
    std::reverse(query.begin(), query.end());
    std::reverse(ref.begin(), ref.end());
    std::reverse(query_ptrs[PTRS].begin(), query_ptrs[PTRS].end());
    std::reverse(query_ptrs[FLAGS].begin(), query_ptrs[FLAGS].end());
    std::reverse(ref_ptrs[PTRS].begin(), ref_ptrs[PTRS].end());
    std::reverse(ref_ptrs[FLAGS].begin(), ref_ptrs[FLAGS].end());

    // update pointers
    for (int & ptr : query_ptrs[PTRS]) { ptr = ref.size()-1 - ptr; }
    for (int & ptr : ref_ptrs[PTRS]) { ptr = query.size()-1 - ptr; }
    for (int & flags : query_ptrs[FLAGS]) { // flip beg/end if only one is set
        if (flags & PTR_VAR_BEG && !(flags & PTR_VAR_END)) {
            flags = PTR_VARIANT | PTR_VAR_END;
        } else if (flags & PTR_VAR_END && !(flags & PTR_VAR_BEG)) {
            flags = PTR_VARIANT | PTR_VAR_BEG;
        }
    }
    for (int & flags : ref_ptrs[FLAGS]) { // flip beg/end if only one is set
        if (flags & PTR_VAR_BEG && !(flags & PTR_VAR_END)) {
            flags = PTR_VARIANT | PTR_VAR_END;
        } else if (flags & PTR_VAR_END && !(flags & PTR_VAR_BEG)) {
            flags = PTR_VARIANT | PTR_VAR_BEG;
        }
    }
    return;
}


/******************************************************************************/


/* Generate the new sequence by applying variants to the reference. */
std::string generate_str(
        std::shared_ptr<fastaData> ref, 
        std::shared_ptr<ctgVariants> vars, const std::string & ctg,
        int var_idx, int end_idx, int beg_pos, int end_pos, int min_qual=0) {

    std::string str = "";
    for (int ref_pos = beg_pos; ref_pos < end_pos; ) {

        // VARIANT
        if (var_idx < end_idx && ref_pos == vars->poss[var_idx]) {
            if (vars->var_quals[var_idx] >= min_qual) {
                switch (vars->types[var_idx]) {
                    case TYPE_INS:
                        str += vars->alts[var_idx];
                        break;
                    case TYPE_DEL:
                        ref_pos += vars->refs[var_idx].size();
                        break;
                    case TYPE_SUB:
                        str += vars->alts[var_idx];
                        ref_pos++;
                        break;
                    case TYPE_CPX:
                        str += vars->alts[var_idx];
                        ref_pos += vars->refs[var_idx].size();
                        break;
                }
            }
            var_idx++; // next variant

        } else { // NO VARIANT
            try {
                // add entire ref seq up to next variant or region end
                int ref_end = (var_idx < end_idx) ? 
                    std::min(end_pos, vars->poss[var_idx]) : end_pos;
                str += ref->fasta.at(ctg).substr(ref_pos, ref_end-ref_pos);
                ref_pos = ref_end;
            } catch (const std::out_of_range & e) {
                ERROR("Contig %s not in reference FASTA (generate_str)", ctg.data());
            }
        }
    }
    return str;
}


/******************************************************************************/


/* For each position on hap1 and hap2, it generates pointers from the query to 
 * the reference and vice versa, as well as the strings for alignment.
 */
void generate_ptrs_strs(
        std::string & query_str, std::string & ref_str,
        std::vector< std::vector<int> > & query_ptrs, 
        std::vector< std::vector<int> > & ref_ptrs,
        std::shared_ptr<ctgVariants> query_vars,
        size_t query_clust_beg_idx, size_t query_clust_end_idx,
        int beg_pos, int end_pos, 
        std::shared_ptr<fastaData> ref, const std::string & ctg
        ) {

    // generate query and ref strings and pointers
    query_ptrs.resize(PTR_DIMS);
    ref_ptrs.resize(PTR_DIMS);
    int query_var_idx = query_vars->clusters.size() ? 
            query_vars->clusters[query_clust_beg_idx] : 0;
    int query_end_idx = query_vars->clusters.size() ? 
            query_vars->clusters[query_clust_end_idx] : 0;
    for (int ref_pos = beg_pos; ref_pos < end_pos; ) {

        if (query_var_idx < query_end_idx && 
                ref_pos == query_vars->poss[query_var_idx]) { // start query variant
            switch (query_vars->types[query_var_idx]) {
                case TYPE_INS:
                    query_ptrs[PTRS].insert(query_ptrs[PTRS].end(), 
                            query_vars->alts[query_var_idx].size(), ref_str.size()-1);
                    query_ptrs[FLAGS].insert(query_ptrs[FLAGS].end(),
                            query_vars->alts[query_var_idx].size(), PTR_VARIANT);
                    query_ptrs[FLAGS][ query_ptrs[FLAGS].size()-1] |= PTR_VAR_END;
                    query_ptrs[FLAGS][ query_ptrs[FLAGS].size() - 
                        query_vars->alts[query_var_idx].size()] |= PTR_VAR_BEG;
                    query_str += query_vars->alts[query_var_idx];
                    break;
                case TYPE_DEL:
                    ref_ptrs[PTRS].insert(ref_ptrs[PTRS].end(),
                            query_vars->refs[query_var_idx].size(), query_str.size()-1);
                    ref_ptrs[FLAGS].insert(ref_ptrs[FLAGS].end(),
                            query_vars->refs[query_var_idx].size(), PTR_VARIANT);
                    ref_ptrs[FLAGS][ ref_ptrs[FLAGS].size()-1] |= PTR_VAR_END;
                    ref_ptrs[FLAGS][ ref_ptrs[FLAGS].size() - 
                        query_vars->refs[query_var_idx].size()] |= PTR_VAR_BEG;
                    ref_str += query_vars->refs[query_var_idx];
                    ref_pos += query_vars->refs[query_var_idx].size();
                    break;
                case TYPE_SUB:
                    ref_ptrs[PTRS].push_back(query_str.size());
                    ref_ptrs[FLAGS].push_back(PTR_VARIANT|PTR_VAR_BEG|PTR_VAR_END);
                    query_ptrs[PTRS].push_back(ref_str.size());
                    query_ptrs[FLAGS].push_back(PTR_VARIANT|PTR_VAR_BEG|PTR_VAR_END);
                    ref_str += query_vars->refs[query_var_idx];
                    query_str += query_vars->alts[query_var_idx];
                    ref_pos++;
                    break;
                default:
                    ERROR("Unexpected variant type '%d' in generate_ptrs_strs()",
                            query_vars->types[query_var_idx]);
            }
            query_var_idx++; // next variant

        } else { // add all matching ref bases

            try {
                // find next position w/o ref match (next var or end)
                int ref_end = (query_var_idx < query_end_idx) ?
                    query_vars->poss[query_var_idx] : end_pos;

                // add pointers, query
                std::vector<int> new_query_ptrs(ref_end - ref_pos, 0);
                query_ptrs[FLAGS].insert(query_ptrs[FLAGS].end(),
                        new_query_ptrs.begin(), new_query_ptrs.end()); // zeros
                for(size_t i = 0; i < new_query_ptrs.size(); i++)
                    new_query_ptrs[i] = ref_str.size() + i;
                query_ptrs[PTRS].insert(query_ptrs[PTRS].end(), 
                        new_query_ptrs.begin(), new_query_ptrs.end());

                // add pointers, ref
                std::vector<int> new_ref_ptrs(ref_end - ref_pos, 0);
                ref_ptrs[FLAGS].insert(ref_ptrs[FLAGS].end(),
                        new_ref_ptrs.begin(), new_ref_ptrs.end()); // zeros
                for(size_t i = 0; i < new_ref_ptrs.size(); i++)
                    new_ref_ptrs[i] = query_str.size() + i;
                ref_ptrs[PTRS].insert(ref_ptrs[PTRS].end(), 
                        new_ref_ptrs.begin(), new_ref_ptrs.end());

                // add sequence, update positions
                std::string matches = ref->fasta.at(ctg).substr(ref_pos, ref_end-ref_pos);
                query_str += matches;
                ref_str += matches;
                ref_pos = ref_end;

            } catch (const std::out_of_range & e) {
                ERROR("Contig '%s' not present in reference FASTA", ctg.data());
            }
        }
    }
}


/******************************************************************************/

/* Calculate initial forward-pass alignment for truth and query strings, given
 * pointers to/from reference. This function generates the pointer matrix, 
 * scores for each alignment, and pointer to if alignment ends on QUERY/REF.
 */
void calc_prec_recall_aln(
        const std::string & query1, const std::string & query2,
        const std::string & truth1, const std::string & truth2, 
        const std::string & ref,
        const std::vector< std::vector<int> > & query1_ref_ptrs, 
        const std::vector< std::vector<int> > & ref_query1_ptrs,
        const std::vector< std::vector<int> > & query2_ref_ptrs, 
        const std::vector< std::vector<int> > & ref_query2_ptrs,
        const std::vector< std::vector<int> > & truth1_ref_ptrs, 
        const std::vector< std::vector<int> > & truth2_ref_ptrs,
        std::vector<int> & s, 
        std::unordered_map<idx1, int> & ptrs,
        std::vector< std::unordered_map<idx1, idx1> > & swap_pred_map,
        std::vector<int> & pr_query_ref_end, bool print
        ) {
    
    // set loop variables
    int ref_len = ref.size();
    std::vector<std::string> query {query1, query1, query2, query2};
    std::vector<std::string> truth {truth1, truth2, truth1, truth2};
    std::vector< std::vector< std::vector<int> > > query_ref_ptrs {
            query1_ref_ptrs, query1_ref_ptrs, query2_ref_ptrs, query2_ref_ptrs };
    std::vector< std::vector< std::vector<int> > > truth_ref_ptrs {
            truth1_ref_ptrs, truth2_ref_ptrs, truth1_ref_ptrs, truth2_ref_ptrs };
    std::vector< std::vector< std::vector<int> > > ref_query_ptrs {
            ref_query1_ptrs, ref_query1_ptrs, ref_query2_ptrs, ref_query2_ptrs };
    std::vector<int> query_lens = 
            {int(query1.size()), int(query1.size()), int(query2.size()), int(query2.size())};
    std::vector<int> truth_lens = 
            {int(truth1.size()), int(truth2.size()), int(truth1.size()), int(truth2.size())};

    std::unordered_set<idx1> done;

    // for each combination of query and truth
    for (int i = 0; i < 4; i++) {
        int qi = 2*i + QUERY; // query index (ptrs)
        int ri = 2*i + REF;   // ref index   (ptrs)
        
        // set first wavefront
        std::queue<idx1> queue; // still to be explored in this wave
        queue.push({qi, 0, 0});
        set(ptrs, idx1(qi, 0, 0), PTR_MAT);
        done.insert(idx1(qi,0,0));
        queue.push({ri, 0, 0});
        set(ptrs, idx1(ri, 0, 0), PTR_MAT);
        done.insert(idx1(ri,0,0));

        // continue looping until full alignment found
        std::unordered_set<idx1> curr_wave; // everything explored this wave
        std::unordered_set<idx1> prev_wave; // everything explored prev wave
        /* if (print) printf("\nFWD %s aln: (%d|%d, %d|%d, %d)\n", aln_strs[i].data(), */ 
        /*         qi, ri, query_lens[i], ref_len, truth_lens[i]); */
        while (true) {
            /* if (print) printf("  s = %d\n", s[i]); */
            if (queue.empty()) ERROR("Empty queue in 'prec_recall_aln()'.");

            // EXTEND WAVEFRONT (stay at same score)
            while (!queue.empty()) {
                idx1 x = queue.front(); queue.pop();
                /* if (print) printf("    x = (%s, %d, %d)\n", */ 
                /*         (x.hi % 2) ? "REF  " : "QUERY", x.qri, x.ti); */
                prev_wave.insert(x);
                if (x.hi == qi) { // QUERY
                    // allow match on query
                    idx1 y(qi, x.qri+1, x.ti+1);
                    if (y.qri < query_lens[i] && y.ti < truth_lens[i] &&
                            query[i][y.qri] == truth[i][y.ti]) {
                        if (!contains(done, y)) {
                            if (!contains(curr_wave, y)) {
                                queue.push(y); curr_wave.insert(y);
                            }
                            set(ptrs, y, PTR_MAT);
                        }
                    }
                    // allow match, swapping to reference
                    idx1 z(ri, query_ref_ptrs[i][PTRS][x.qri]+1, x.ti+1);
                    if ( (!(query_ref_ptrs[i][FLAGS][x.qri] & PTR_VARIANT) ||
                            query_ref_ptrs[i][FLAGS][x.qri] & PTR_VAR_END) &&
                         (!(truth_ref_ptrs[i][FLAGS][x.ti] & PTR_VARIANT) ||
                            truth_ref_ptrs[i][FLAGS][x.ti] & PTR_VAR_END)) {
                        if (z.qri < ref_len && z.ti < truth_lens[i] &&
                                ref[z.qri] == truth[i][z.ti]) {
                            if (!contains(done, z)) {
                                if (!contains(curr_wave, z)) {
                                    queue.push(z); curr_wave.insert(z);
                                }
                                set(ptrs, z, PTR_SWP_MAT);
                                swap_pred_map[i][z] = x;
                            }
                        }
                    }
                } else { // x.hi == ri == REF
                    // allow match
                    idx1 y(ri, x.qri+1, x.ti+1);
                    if (y.qri < ref_len && y.ti < truth_lens[i] &&
                            ref[y.qri] == truth[i][y.ti]) {
                        if (!contains(done, y)) {
                            if (!contains(curr_wave, y)) {
                                queue.push(y); curr_wave.insert(y);
                            }
                            set(ptrs, y, PTR_MAT);
                        }
                    }
                    // allow match, swapping to query
                    idx1 z(qi, ref_query_ptrs[i][PTRS][x.qri]+1, x.ti+1);
                    if ( (!(ref_query_ptrs[i][FLAGS][x.qri] & PTR_VARIANT) ||
                            ref_query_ptrs[i][FLAGS][x.qri] & PTR_VAR_END) &&
                         (!(truth_ref_ptrs[i][FLAGS][x.ti] & PTR_VARIANT) ||
                            truth_ref_ptrs[i][FLAGS][x.ti] & PTR_VAR_END)) {
                        if (z.qri < query_lens[i] && z.ti < truth_lens[i] &&
                                query[i][z.qri] == truth[i][z.ti]) {
                            if (!contains(done, z)) {
                                if (!contains(curr_wave, z)) {
                                    queue.push(z); curr_wave.insert(z);
                                }
                                set(ptrs, z, PTR_SWP_MAT);
                                swap_pred_map[i][z] = x;
                            }
                        }
                    }
                }
            }

            // mark all cells visited this wave as done
            for (auto x : curr_wave) { done.insert(x); }
            curr_wave.clear();

            // exit if we're done aligning
            if (contains(done, idx1(qi, query_lens[i]-1, truth_lens[i]-1)) ||
                contains(done, idx1(ri, ref_len-1, truth_lens[i]-1))) break;


            // NEXT WAVEFRONT (increase score by one)
            for (auto x : prev_wave) {
                int qr_len = (x.hi == qi) ? query_lens[i] : ref_len;
                if (x.qri+1 < qr_len) { // INS
                    idx1 y(x.hi, x.qri+1, x.ti);
                    if (!contains(done, y) && !contains(curr_wave, y)) {
                        queue.push(y);
                        curr_wave.insert(y);
                    }
                    if (!contains(done, y))
                        set(ptrs, y, PTR_INS);
                }
                if (x.ti+1 < truth_lens[i]) { // DEL
                    idx1 y(x.hi, x.qri, x.ti+1);
                    if (!contains(done, y) && !contains(curr_wave, y)) {
                        queue.push(y);
                        curr_wave.insert(y);
                    }
                    if (!contains(done, y))
                        set(ptrs, y, PTR_DEL);
                }
                if (x.qri+1 < qr_len && x.ti+1 < truth_lens[i]) { // SUB
                    idx1 y(x.hi, x.qri+1, x.ti+1);
                    if (!contains(done, y) && !contains(curr_wave, y)) {
                        queue.push(y);
                        curr_wave.insert(y);
                    }
                    if (!contains(done, y))
                        set(ptrs, y, PTR_SUB);
                }
            }
            prev_wave.clear();
            s[i]++;
        } // while loop (this alignment)

        if (print) printf("\nAlignment %s, aln_ptrs\n", aln_strs[i].data());
        if (print) printf("\nQUERY");
        if (print) print_ptrs(ptrs, qi, query[i], truth[i]);
        if (print) printf("\nREF");
        if (print) print_ptrs(ptrs, ri, ref, truth[i]);

        // save where to start backtrack (prefer ref: omit vars which don't reduce ED)
        if (contains(done, idx1(ri, ref_len-1, truth_lens[i]-1))) {
            pr_query_ref_end[i] = ri;
        } else if (contains(done, idx1(qi, query_lens[i]-1, truth_lens[i]-1))) {
            pr_query_ref_end[i] = qi;
        } else { ERROR("Alignment not finished in 'prec_recall_aln()'."); }

    } // 4 alignments
} // function


/******************************************************************************/


int store_phase( 
        std::shared_ptr<superclusterData> clusterdata_ptr, 
        const std::string & ctg,
        const std::vector<int> & s
        ) {

    // calculate best phasing
    int orig_phase_dist = s[QUERY1_TRUTH1] + s[QUERY2_TRUTH2];
    int swap_phase_dist = s[QUERY2_TRUTH1] + s[QUERY1_TRUTH2];
    int phase = PHASE_NONE; // default either way if equal dist
    if (orig_phase_dist < swap_phase_dist) phase = PHASE_ORIG;
    if (swap_phase_dist < orig_phase_dist) phase = PHASE_SWAP;

    // save alignment information
    clusterdata_ptr->ctg_superclusters[ctg]->add_phasing(
            phase, orig_phase_dist, swap_phase_dist);
    return phase;
}


/******************************************************************************/


/* Given alignment pointers (aln_ptrs) from minED forward-pass, as well as 
 * reference pointers to/from truth/query, backtrack and find path which 
 * maximizes FP variants. Results are stored in path_ptrs. Then call 
 * get_prec_recall_path_sync() to parse this path.
 */
void calc_prec_recall_path(
        const std::string & ref,
        const std::string & query1, const std::string & query2, 
        const std::string & truth1, const std::string & truth2, 
        std::vector< std::vector<idx1> > & path, 
        std::vector< std::vector<bool> > & sync, 
        std::vector< std::vector<bool> > & edits, 
        std::unordered_map<idx1, int> & aln_ptrs,
        std::unordered_map<idx1, int> & path_ptrs,
        std::unordered_map<idx1, int> & path_scores,
        const std::vector< std::vector<int> > & query1_ref_ptrs, 
        const std::vector< std::vector<int> > & ref_query1_ptrs,
        const std::vector< std::vector<int> > & query2_ref_ptrs, 
        const std::vector< std::vector<int> > & ref_query2_ptrs,
        const std::vector< std::vector<int> > & truth1_ref_ptrs, 
        const std::vector< std::vector<int> > & truth2_ref_ptrs,
        const std::vector< std::unordered_map<idx1,idx1> > & swap_pred_map,
        const std::vector<int> & pr_query_ref_end, bool print
        ) {

    std::vector<std::string> query = {query1, query1, query2, query2};
    std::vector<std::string> truth = {truth1, truth2, truth1, truth2};
    std::vector< std::vector< std::vector<int> > > query_ref_ptrs = { 
            query1_ref_ptrs, query1_ref_ptrs, query2_ref_ptrs, query2_ref_ptrs };
    std::vector< std::vector< std::vector<int> > > ref_query_ptrs = { 
            ref_query1_ptrs, ref_query1_ptrs, ref_query2_ptrs, ref_query2_ptrs };
    std::vector< std::vector< std::vector<int> > > truth_ref_ptrs = { 
            truth1_ref_ptrs, truth2_ref_ptrs, truth1_ref_ptrs, truth2_ref_ptrs };
    std::vector<int> pr_query_ref_beg(4);
    std::vector< std::vector<bool> > ref_loc_sync;
    std::unordered_set<idx1> done;

    for (int i = 0; i < 4; i++) {
        /* printf("\nBKWD %s path\n", aln_strs[i].data()); */

        // init
        int qi = i*2 + QUERY;
        int ri = i*2 + REF;
        path.push_back(std::vector<idx1>());
        sync.push_back(std::vector<bool>());
        edits.push_back(std::vector<bool>());
        ref_loc_sync.push_back(std::vector<bool>(ref_query_ptrs[i][0].size(), true));

        // backtrack start
        std::queue<idx1> queue;
        int start_hi = pr_query_ref_end[i];
        int start_qri = (start_hi % 2 == QUERY) ? query_ref_ptrs[i][0].size()-1 :
                ref_query_ptrs[i][0].size()-1;
        int start_ti = truth_ref_ptrs[i][0].size()-1;
        idx1 start(start_hi, start_qri, start_ti);
        path_ptrs[start] = PTR_MAT;
        set(aln_ptrs, start, MAIN_PATH); // all alignments end here
        path_scores[start] = 0;
        queue.push(start);
        std::unordered_set<idx1> curr_wave, prev_wave;

        while (true) {
            while (!queue.empty()) {
                idx1 x = queue.front(); queue.pop(); // current cell
                prev_wave.insert(x);
                /* printf("(%d, %d, %d)\n", x.hi, x.qri, x.ti); */

                // MATCH movement
                if (test(aln_ptrs, x, PTR_MAT) && x.qri > 0 && x.ti > 0) {

                    // get next cell, add to path
                    idx1 y = idx1(x.hi, x.qri-1, x.ti-1);
                    set(aln_ptrs, y, PATH);

                    // check if mvmt is in variant
                    bool in_truth_var = truth_ref_ptrs[i][FLAGS][x.ti] & PTR_VARIANT;
                    in_truth_var &= !(truth_ref_ptrs[i][FLAGS][x.ti] & PTR_VAR_BEG);
                    bool in_query_var = (x.hi == ri) ? false : 
                        query_ref_ptrs[i][FLAGS][x.qri] & PTR_VARIANT;
                    in_query_var &= !(query_ref_ptrs[i][FLAGS][x.qri] & PTR_VAR_BEG);

                    // if off-diagonal or in truth var, this ref loc isn't a sync point
                    if (y.hi == ri && (
                            in_truth_var || y.qri != truth_ref_ptrs[i][PTRS][y.ti]))
                        ref_loc_sync[i][y.qri] = false;
                    if (y.hi == qi && !in_query_var && ( in_truth_var || 
                            query_ref_ptrs[i][PTRS][y.qri] != truth_ref_ptrs[i][PTRS][y.ti]))
                        ref_loc_sync[i][ query_ref_ptrs[i][PTRS][y.qri] ] = false;

                    // check for fp
                    bool is_fp = x.hi == ri && ((ref_query_ptrs[i][PTRS][x.qri] != 
                            ref_query_ptrs[i][PTRS][y.qri]+1) || // INS
                            ref_query_ptrs[i][FLAGS][x.qri] & PTR_VAR_BEG); // SUB/DEL

                    // update score
                    if (!contains(done, y) &&
                            (!contains(path_scores, y) ||
                            path_scores[x] + is_fp > path_scores[y])) {
                        path_ptrs[y] = PTR_MAT;
                        path_scores[y] = path_scores[x] + is_fp;
                        if (!contains(curr_wave, y)) curr_wave.insert(y);
                        queue.push(y);
                    } else if (!contains(done, y) &&
                            path_scores[x] + is_fp == path_scores[y]) {
                        set(path_ptrs, y, PTR_MAT);
                    }

                    if (print) printf("(%s, %2d, %2d) --%s-> (%s, %2d, %2d) %s %s %s %s\n",
                            x.hi % 2 ? "REF" : "QRY", x.qri, x.ti, "MAT",
                            y.hi % 2 ? "REF" : "QRY", y.qri, y.ti,
                            in_query_var ? GREEN("QUERY_VAR").data() : RED("QUERY_VAR").data(), 
                            in_truth_var ? GREEN("TRUTH_VAR").data() : RED("TRUTH_VAR").data(),
                            ref_loc_sync[i][y.hi == ri ? y.qri : query_ref_ptrs[i][PTRS][y.qri]] ? GREEN("REF_SYNC").data() : RED("REF_SYNC").data(),
                            is_fp ? GREEN("FP").data() : RED("FP").data());

                }


                if (x.hi == ri) { // REF -> QUERY SWAP mvmt
                    if (test(aln_ptrs, x, PTR_SWP_MAT) &&
                            (!(ref_query_ptrs[i][FLAGS][x.qri] & PTR_VARIANT) ||
                            ref_query_ptrs[i][FLAGS][x.qri] & PTR_VAR_BEG) && 
                            x.qri > 0 && x.ti > 0) {

                        // get next cell, add to path
                        if (!contains(swap_pred_map[i], x))
                            ERROR("No swap predecessor, but PTR_SWP_MAT set.");
                        idx1 z = swap_pred_map[i].at(x);
                        set(aln_ptrs, z, PATH);

                        // check if mvmt is in variant
                        bool in_truth_var = truth_ref_ptrs[i][FLAGS][x.ti] & PTR_VARIANT;
                        in_truth_var &= !(truth_ref_ptrs[i][FLAGS][x.ti] & PTR_VAR_BEG);
                        bool in_query_var = (x.hi == ri) ? false : 
                            query_ref_ptrs[i][FLAGS][x.qri] & PTR_VARIANT;
                        in_query_var &= !(query_ref_ptrs[i][FLAGS][x.qri] & PTR_VAR_BEG);

                        // if off-diag or in truth var, this ref loc isn't a sync point
                        if (!in_query_var && ( in_truth_var || 
                                query_ref_ptrs[i][PTRS][z.qri] != truth_ref_ptrs[i][PTRS][z.ti]))
                            ref_loc_sync[i][query_ref_ptrs[i][PTRS][z.qri]] = false;

                        // check for fp
                        bool is_fp = x.hi == ri && ((ref_query_ptrs[i][PTRS][x.qri] != 
                                ref_query_ptrs[i][PTRS][x.qri-1]+1) || // INS
                                ref_query_ptrs[i][FLAGS][x.qri] & PTR_VAR_BEG); // SUB/DEL

                        // update score
                        if (!contains(done, z) &&
                            (!contains(path_scores, z) ||
                                path_scores[x] + is_fp > path_scores[z])) {
                            path_ptrs[z] = PTR_SWP_MAT;
                            path_scores[z] = path_scores[x] + is_fp;
                            if (!contains(curr_wave, z)) curr_wave.insert(z);
                            queue.push(z);
                        } else if (!contains(done, z) &&
                                path_scores[x] + is_fp == path_scores[z]) {
                            set(path_ptrs, z, PTR_SWP_MAT);
                        }

                        if (print) printf("(%s, %2d, %2d) --%s-> (%s, %2d, %2d) %s %s %s %s\n",
                                x.hi % 2 ? "REF" : "QRY", x.qri, x.ti, "MAT",
                                z.hi % 2 ? "REF" : "QRY", z.qri, z.ti,
                                in_query_var ? GREEN("QUERY_VAR").data() : RED("QUERY_VAR").data(), 
                                in_truth_var ? GREEN("TRUTH_VAR").data() : RED("TRUTH_VAR").data(),
                                ref_loc_sync[i][z.hi == ri ? z.qri : query_ref_ptrs[i][PTRS][z.qri]] ? GREEN("REF_SYNC").data() : RED("REF_SYNC").data(),
                                is_fp ? GREEN("FP").data() : RED("FP").data());
                    }

                } else { // QUERY -> REF SWAP mvmt
                    if (test(aln_ptrs, x, PTR_SWP_MAT) &&
                            (!(query_ref_ptrs[i][FLAGS][x.qri] & PTR_VARIANT) ||
                            query_ref_ptrs[i][FLAGS][x.qri] & PTR_VAR_BEG) && 
                            x.qri > 0 && x.ti > 0) {

                        // add to path
                        if (!contains(swap_pred_map[i], x))
                            ERROR("No swap predecessor, but PTR_SWP_MAT set.");
                        idx1 z = swap_pred_map[i].at(x);
                        set(aln_ptrs, z, PATH);

                        // check if mvmt is in variant
                        bool in_truth_var = truth_ref_ptrs[i][FLAGS][x.ti] & PTR_VARIANT;
                        in_truth_var &= !(truth_ref_ptrs[i][FLAGS][x.ti] & PTR_VAR_BEG);
                        bool in_query_var = (x.hi == ri) ? false : 
                            query_ref_ptrs[i][FLAGS][x.qri] & PTR_VARIANT;
                        in_query_var &= !(query_ref_ptrs[i][FLAGS][x.qri] & PTR_VAR_BEG);

                        // if off-diag or in truth var, this ref loc isn't a sync point
                        if (in_truth_var || z.qri != truth_ref_ptrs[i][PTRS][z.ti])
                            ref_loc_sync[i][z.qri] = false;

                        // no need to check for FP since we were on QUERY, not REF

                        // update score
                        if (!contains(done, z) &&
                            (!contains(path_scores, z) ||
                                path_scores[x] > path_scores[z])) {
                            path_ptrs  [z] = PTR_SWP_MAT;
                            path_scores[z] = path_scores[x];
                            if (!contains(curr_wave, z)) curr_wave.insert(z);
                            queue.push(z);
                        } else if (!contains(done, z) &&
                                path_scores[x] == path_scores[z]) {
                            set(path_ptrs, z, PTR_SWP_MAT);
                        }

                        if (print) printf("(%s, %2d, %2d) --%s-> (%s, %2d, %2d) %s %s %s\n",
                                x.hi % 2 ? "REF" : "QRY", x.qri, x.ti, "MAT",
                                z.hi % 2 ? "REF" : "QRY", z.qri, z.ti,
                                in_query_var ? GREEN("QUERY_VAR").data() : RED("QUERY_VAR").data(), 
                                in_truth_var ? GREEN("TRUTH_VAR").data() : RED("TRUTH_VAR").data(),
                                ref_loc_sync[i][z.hi == ri ? z.qri : query_ref_ptrs[i][PTRS][z.qri]] ? GREEN("REF_SYNC").data() : RED("REF_SYNC").data());
                    }
                }

            } // end of curr_wave, same score

            for (idx1 x : curr_wave) done.insert(x);
            curr_wave.clear();
            if (contains(done, idx1(qi, 0, 0)) || contains(done, idx1(ri, 0, 0))) break;

            for (idx1 x : prev_wave) {

                // SUB movement
                if (test(aln_ptrs, x, PTR_SUB) &&
                        x.qri > 0 && x.ti > 0) {

                    // get next cell, add to path
                    idx1 y = idx1(x.hi, x.qri-1, x.ti-1);
                    set(aln_ptrs, y, PATH);

                    // check if mvmt is in variant
                    bool in_truth_var = truth_ref_ptrs[i][FLAGS][x.ti] & PTR_VARIANT;
                    in_truth_var &= !(truth_ref_ptrs[i][FLAGS][x.ti] & PTR_VAR_BEG);
                    bool in_query_var = (x.hi == ri) ? false : 
                        query_ref_ptrs[i][FLAGS][x.qri] & PTR_VARIANT;
                    in_query_var &= !(query_ref_ptrs[i][FLAGS][x.qri] & PTR_VAR_BEG);

                    // if off-diagonal or in truth var, this ref loc isn't a sync point
                    if (y.hi == ri && (
                            in_truth_var || y.qri != truth_ref_ptrs[i][PTRS][y.ti]))
                        ref_loc_sync[i][y.qri] = false;
                    if (y.hi == qi && !in_query_var && ( in_truth_var || 
                            query_ref_ptrs[i][PTRS][y.qri] != truth_ref_ptrs[i][PTRS][y.ti]))
                        ref_loc_sync[i][query_ref_ptrs[i][PTRS][y.qri]] = false;

                    // check for fp
                    bool is_fp = x.hi == ri && ((ref_query_ptrs[i][PTRS][x.qri] != 
                            ref_query_ptrs[i][PTRS][y.qri]+1) || // INS
                            ref_query_ptrs[i][FLAGS][x.qri] & PTR_VAR_BEG); // SUB/DEL

                    // update score
                    if (!contains(done, y) &&
                            (!contains(path_scores, y) ||
                            path_scores[x] + is_fp > path_scores[y])) {
                        path_ptrs  [y] = PTR_SUB;
                        path_scores[y] = path_scores[x] + is_fp;
                        if (!contains(curr_wave, y)) curr_wave.insert(y);
                        queue.push(y);
                    } else if (!contains(done, y) &&
                            path_scores[x] + is_fp == path_scores[y]) {
                        set(path_ptrs, y, PTR_SUB);
                    }

                    if (print) printf("(%s, %2d, %2d) --%s-> (%s, %2d, %2d) %s %s %s %s\n",
                            x.hi % 2 ? "REF" : "QRY", x.qri, x.ti, "MAT",
                            y.hi % 2 ? "REF" : "QRY", y.qri, y.ti,
                            in_query_var ? GREEN("QUERY_VAR").data() : RED("QUERY_VAR").data(), 
                            in_truth_var ? GREEN("TRUTH_VAR").data() : RED("TRUTH_VAR").data(),
                            ref_loc_sync[i][y.hi == ri ? y.qri : query_ref_ptrs[i][PTRS][y.qri]] ? GREEN("REF_SYNC").data() : RED("REF_SYNC").data(),
                            is_fp ? GREEN("FP").data() : RED("FP").data());
                }

                // INS movement
                if (test(aln_ptrs, x, PTR_INS) && x.qri > 0) {

                    // get next cell, add to path
                    idx1 y = idx1(x.hi, x.qri-1, x.ti);
                    set(aln_ptrs, y, PATH);

                    // check if mvmt is in variant
                    bool in_truth_var = truth_ref_ptrs[i][FLAGS][x.ti] & PTR_VARIANT;
                    bool in_query_var = (x.hi == ri) ? false : 
                        query_ref_ptrs[i][FLAGS][x.qri] & PTR_VARIANT;
                    in_query_var &= !(query_ref_ptrs[i][FLAGS][x.qri] & PTR_VAR_BEG);

                    // if off-diagonal or in truth var, this ref loc isn't a sync point
                    if (y.hi == ri && (
                            in_truth_var || y.qri != truth_ref_ptrs[i][PTRS][y.ti]))
                        ref_loc_sync[i][y.qri] = false;
                    if (y.hi == qi && !in_query_var && ( in_truth_var || 
                            query_ref_ptrs[i][PTRS][y.qri] != truth_ref_ptrs[i][PTRS][y.ti]))
                        ref_loc_sync[i][query_ref_ptrs[i][PTRS][y.qri]] = false;

                    // check for fp
                    bool is_fp = x.hi == ri && ((ref_query_ptrs[i][PTRS][x.qri] != 
                            ref_query_ptrs[i][PTRS][y.qri]+1) || // INS
                            ref_query_ptrs[i][FLAGS][x.qri] & PTR_VAR_BEG); // SUB/DEL

                    // update score
                    if (!contains(done, y) &&
                            (!contains(path_scores, y) ||
                            path_scores[x] + is_fp > path_scores[y])) {
                        path_ptrs  [y] = PTR_INS;
                        path_scores[y] = path_scores[x] + is_fp;
                        if (!contains(curr_wave, y)) curr_wave.insert(y);
                        queue.push(y);
                    } else if (!contains(done, y) &&
                            path_scores[x] + is_fp == path_scores[y]) {
                        set(path_ptrs, y, PTR_INS);
                    }

                    if (print) printf("(%s, %2d, %2d) --%s-> (%s, %2d, %2d) %s %s %s %s\n",
                            x.hi % 2 ? "REF" : "QRY", x.qri, x.ti, "MAT",
                            y.hi % 2 ? "REF" : "QRY", y.qri, y.ti,
                            in_query_var ? GREEN("QUERY_VAR").data() : RED("QUERY_VAR").data(), 
                            in_truth_var ? GREEN("TRUTH_VAR").data() : RED("TRUTH_VAR").data(),
                            ref_loc_sync[i][y.hi == ri ? y.qri : query_ref_ptrs[i][PTRS][y.qri]] ? GREEN("REF_SYNC").data() : RED("REF_SYNC").data(),
                            is_fp ? GREEN("FP").data() : RED("FP").data());
                }

                // DEL movement
                if (test(aln_ptrs, x, PTR_DEL) && x.ti > 0) {

                    // add to path
                    idx1 y = idx1(x.hi, x.qri, x.ti-1);
                    set(aln_ptrs, y, PATH);

                    bool in_truth_var = truth_ref_ptrs[i][FLAGS][x.ti] & PTR_VARIANT;
                    in_truth_var &= !(truth_ref_ptrs[i][FLAGS][x.ti] & PTR_VAR_BEG);
                    bool in_query_var = (x.hi == ri) ? false : 
                        query_ref_ptrs[i][FLAGS][x.qri] & PTR_VARIANT;

                    // no need to check for FP since we don't consume a REF/QUERY base
                    
                    // update score
                    if (!contains(done, y) &&
                            (!contains(path_scores, y) ||
                            path_scores[x] > path_scores[y])) {
                        path_ptrs  [y] = PTR_DEL;
                        path_scores[y] = path_scores[x];
                        if (!contains(curr_wave, y)) curr_wave.insert(y);
                        queue.push(y);
                    } else if (!contains(done, y) &&
                            path_scores[x] == path_scores[y]) {
                        set(path_ptrs, y, PTR_DEL);
                    }

                    if (print) printf("(%s, %2d, %2d) --%s-> (%s, %2d, %2d) %s %s %s\n",
                            x.hi % 2 ? "REF" : "QRY", x.qri, x.ti, "MAT",
                            y.hi % 2 ? "REF" : "QRY", y.qri, y.ti,
                            in_query_var ? GREEN("QUERY_VAR").data() : RED("QUERY_VAR").data(), 
                            in_truth_var ? GREEN("TRUTH_VAR").data() : RED("TRUTH_VAR").data(),
                            ref_loc_sync[i][y.hi == ri ? y.qri : query_ref_ptrs[i][PTRS][y.qri]] ? GREEN("REF_SYNC").data() : RED("REF_SYNC").data());
                }

            } // done adding to next wave
            prev_wave.clear();
        }

        // set pointers
        if (test(aln_ptrs, idx1(ri, 0, 0), PATH))
            pr_query_ref_beg[i] = ri;
        else
            pr_query_ref_beg[i] = qi;

        // debug print
        if (print) printf("Alignment %s, path_ptrs\n", aln_strs[i].data());
        if (print) printf("\nQUERY");
        if (print) print_ptrs(path_ptrs, qi, query[i], truth[i]);
        if (print) printf("\nREF");
        if (print) print_ptrs(path_ptrs, ri, ref, truth[i]);

    } // 4 alignments

    // get path and sync points
    get_prec_recall_path_sync(path, sync, edits, aln_ptrs, path_ptrs, ref_loc_sync,
            query1_ref_ptrs, ref_query1_ptrs, query2_ref_ptrs, ref_query2_ptrs,
            truth1_ref_ptrs, truth2_ref_ptrs, pr_query_ref_beg, print
    );

}


/******************************************************************************/


/* Follow path which maximizes FP and minimizes ED, saving sync points and edits.
 */
void get_prec_recall_path_sync(
        std::vector< std::vector<idx1> > & path, 
        std::vector< std::vector<bool> > & sync, 
        std::vector< std::vector<bool> > & edits, 
        std::unordered_map<idx1, int> & aln_ptrs, 
        std::unordered_map<idx1, int> & path_ptrs, 
        const std::vector< std::vector<bool> > & ref_loc_sync, 
        const std::vector< std::vector<int> > & query1_ref_ptrs, 
        const std::vector< std::vector<int> > & ref_query1_ptrs,
        const std::vector< std::vector<int> > & query2_ref_ptrs, 
        const std::vector< std::vector<int> > & ref_query2_ptrs,
        const std::vector< std::vector<int> > & truth1_ref_ptrs, 
        const std::vector< std::vector<int> > & truth2_ref_ptrs,
        const std::vector<int> & pr_query_ref_beg, bool print
        ) {

    // query <-> ref pointers
    std::vector< std::vector< std::vector<int> > > query_ref_ptrs = { 
            query1_ref_ptrs, query1_ref_ptrs, query2_ref_ptrs, query2_ref_ptrs };
    std::vector< std::vector< std::vector<int> > > truth_ref_ptrs = { 
            truth1_ref_ptrs, truth2_ref_ptrs, truth1_ref_ptrs, truth2_ref_ptrs };
    std::vector< std::vector< std::vector<int> > > ref_query_ptrs = { 
            ref_query1_ptrs, ref_query1_ptrs, ref_query2_ptrs, ref_query2_ptrs };

    for (int i = 0; i < 4; i++) {

        if (print) printf("Alignment %s:\n", aln_strs[i].data());

        if (print) printf("REF LOC SYNC:");
        for(int j = 0; j < int(ref_loc_sync[i].size()); j++) {
            if (print) printf(ref_loc_sync[i][j] ? "=" : "X");
        }
        if (print) printf("\n");

        // init
        int qi = i*2 + QUERY;
        int ri = i*2 + REF;

        // path start
        int hi = pr_query_ref_beg[i];
        int qri = 0;
        int ti = 0;
        int prev_hi, prev_qri, prev_ti = 0;
        int ptr_type = PTR_NONE;

        int r_size = ref_query_ptrs[i][PTRS].size();
        int q_size = query_ref_ptrs[i][PTRS].size();
        int t_size = truth_ref_ptrs[i][PTRS].size();

        // first position is sync point
        path[i].push_back(idx1(hi, qri, ti));
        set(aln_ptrs, idx1(hi, qri, ti), MAIN_PATH);
        set(aln_ptrs, idx1(hi, qri, ti), PTR_SYNC);
        set(path_ptrs, idx1(hi, qri, ti), MAIN_PATH);

        // follow best-path pointers
        while ( (hi == ri && qri < r_size-1) || (hi == qi && qri < q_size-1) || ti < t_size-1) {
            prev_hi = hi; prev_qri = qri; prev_ti = ti;
            if (hi == qi && test(path_ptrs, idx1(hi, qri, ti), PTR_SWP_MAT)) { // prefer FPs
                ptr_type = PTR_SWP_MAT;
                hi = ri;
                qri = query_ref_ptrs[i][PTRS][qri];
                qri++; ti++; edits[i].push_back(false);

            } else if (test(path_ptrs, idx1(hi, qri, ti), PTR_MAT)) {
                ptr_type = PTR_MAT;
                qri++; ti++; edits[i].push_back(false);

            } else if (test(path_ptrs, idx1(hi, qri, ti), PTR_SUB)) {
                ptr_type = PTR_SUB;
                qri++; ti++; edits[i].push_back(true);

            } else if (test(path_ptrs, idx1(hi, qri, ti), PTR_INS)) {
                ptr_type = PTR_INS;
                qri++; edits[i].push_back(true);

            } else if (test(path_ptrs, idx1(hi, qri, ti), PTR_DEL)) {
                ptr_type = PTR_DEL;
                ti++; edits[i].push_back(true);

            } else if (hi == ri && test(path_ptrs, idx1(hi, qri, ti), PTR_SWP_MAT)) {
                ptr_type = PTR_SWP_MAT;
                hi = qi;
                qri = ref_query_ptrs[i][PTRS][qri];
                qri++; ti++; edits[i].push_back(false);

            } else {
                ERROR("No valid pointer (value '%d') at (%d, %d, %d)", 
                        path_ptrs[idx1(hi, qri, ti)], hi, qri, ti);
            }

            if ((hi == qi && qri >= q_size) || (hi == ri && qri >= r_size) || ti >= t_size) break;

            // add point to path
            path[i].push_back(idx1(hi, qri, ti));
            set(aln_ptrs, idx1(hi, qri, ti), MAIN_PATH);
            set(path_ptrs, idx1(hi, qri, ti), MAIN_PATH);

            // set boolean flags for if in query/truth variants
            bool in_truth_var = truth_ref_ptrs[i][FLAGS][ti] & PTR_VARIANT;
            bool in_query_var = (hi == ri) ? false : 
                    query_ref_ptrs[i][FLAGS][qri] & PTR_VARIANT;
            switch(ptr_type) { // not really in variant if just moved into it
                case PTR_MAT: // just consumed truth and query
                case PTR_SWP_MAT:
                case PTR_SUB:
                    in_truth_var &= !(truth_ref_ptrs[i][FLAGS][ti] & PTR_VAR_BEG);
                    if (in_query_var) in_query_var &= !(query_ref_ptrs[i][FLAGS][qri] & PTR_VAR_BEG);
                    break;
                case PTR_INS: // just consumed query
                    if (in_query_var) in_query_var &= !(query_ref_ptrs[i][FLAGS][qri] & PTR_VAR_BEG);
                    break;
                case PTR_DEL: // just consumed truth
                    in_truth_var &= !(truth_ref_ptrs[i][FLAGS][ti] & PTR_VAR_BEG);
                    break;
                default:
                    ERROR("Unexpected pointer type (%d) in calc_path()", ptr_type);
            }

            // SYNC POINTS
            bool prev_sync;
            // no info about movements from previous rows, default to false
            if (prev_qri == 0 || (prev_hi == qi && query_ref_ptrs[i][PTRS][prev_qri] == 0)) {
                prev_sync = false;
            } else { // use ref_loc_sync info, enforce diag, not in variant
                prev_sync = !in_truth_var && !in_query_var && 
                    truth_ref_ptrs[i][PTRS][prev_ti] == (prev_hi == ri ? prev_qri :
                            query_ref_ptrs[i][PTRS][prev_qri]) &&
                    ptr_type & (PTR_MAT | PTR_SWP_MAT | PTR_SUB) &&
                    ref_loc_sync[i][ prev_hi == ri ? prev_qri-1 : 
                            query_ref_ptrs[i][PTRS][prev_qri]-1 ];
            }
            if (prev_qri == 0 && prev_ti == 0) prev_sync = true; // only start
            if (prev_sync) set(aln_ptrs, idx1(prev_hi, prev_qri, prev_ti), PTR_SYNC);
            sync[i].push_back(prev_sync);

            // debug print
            std::string ptr_str;
            if (ptr_type == PTR_MAT) ptr_str = "MAT";
            if (ptr_type == PTR_SUB) ptr_str = "SUB";
            if (ptr_type == PTR_INS) ptr_str = "INS";
            if (ptr_type == PTR_DEL) ptr_str = "DEL";
            if (ptr_type == PTR_SWP_MAT) ptr_str = "SWP";
            if (print) printf("%s (%s, %2d, %2d) --%s-> (%s, %2d, %2d)\n",
                    prev_sync ? "SYNC" : "    ", prev_hi % 2 ? "REF" : "QRY", 
                    prev_qri, prev_ti, 
                    ptr_str.data(), hi % 2 ? "REF" : "QRY", qri, ti);
        }

        // last position is sync
        if (print) printf("SYNC (%s, %2d, %2d) DONE \n",
                hi % 2 ? "REF" : "QRY", qri, ti);
        sync[i].push_back(true);
        set(aln_ptrs, idx1(hi, qri, ti), PTR_SYNC);
    }
}


/******************************************************************************/


void calc_prec_recall(
        std::shared_ptr<superclusterData> clusterdata_ptr, int sc_idx, 
        const std::string & ctg, 
        const std::string & ref,
        const std::string & query1, const std::string & query2, 
        const std::string & truth1, const std::string & truth2, 
        const std::vector< std::vector<idx1> > & path,
        const std::vector< std::vector<bool> > & sync,
        const std::vector< std::vector<bool> > & edits,
        const std::unordered_map<idx1, int> & aln_ptrs, 
        const std::unordered_map<idx1, int> & path_ptrs, 
        const std::vector< std::vector<int> > & query1_ref_ptrs, 
        const std::vector< std::vector<int> > & ref_query1_ptrs,
        const std::vector< std::vector<int> > & query2_ref_ptrs, 
        const std::vector< std::vector<int> > & ref_query2_ptrs,
        const std::vector< std::vector<int> > & truth1_ref_ptrs, 
        const std::vector< std::vector<int> > & truth2_ref_ptrs,
        const std::vector<int> & pr_query_ref_end, int phase, 
        bool print
        ) {

    // set query/truth strings and pointers
    int beg = clusterdata_ptr->ctg_superclusters[ctg]->begs[sc_idx];
    std::vector<std::string> query = {query1, query1, query2, query2};
    std::vector<std::string> truth = {truth1, truth2, truth1, truth2};
    std::vector< std::vector< std::vector<int> > > query_ref_ptrs = { 
            query1_ref_ptrs, query1_ref_ptrs, query2_ref_ptrs, query2_ref_ptrs };
    std::vector< std::vector< std::vector<int> > > ref_query_ptrs = { 
            ref_query1_ptrs, ref_query1_ptrs, ref_query2_ptrs, ref_query2_ptrs };
    std::vector< std::vector< std::vector<int> > > truth_ref_ptrs = { 
            truth1_ref_ptrs, truth2_ref_ptrs, truth1_ref_ptrs, truth2_ref_ptrs };

    // indices into ptr/off matrices depend on decided phasing
    std::vector<int> indices;
    if (phase == PHASE_SWAP) {
        indices.push_back(QUERY1_TRUTH2);
        indices.push_back(QUERY2_TRUTH1);
    } else if (phase == PHASE_ORIG || phase == PHASE_NONE) { // keep
        indices.push_back(QUERY1_TRUTH1);
        indices.push_back(QUERY2_TRUTH2);
    } else {
        ERROR("Unexpected phase (%d)", phase);
    }

    // for only the selected phasing
    for (int i : indices) {

        int ri = i*2 + REF;   // ref index
        int qi = i*2 + QUERY; // call index
        int qhi = i >> 1;     // query hap index
        int thi = i & 1;      // truth hap index

        // set variant ranges
        std::shared_ptr<ctgVariants> query_vars = 
                clusterdata_ptr->ctg_superclusters[ctg]->ctg_variants[QUERY][qhi];
        std::shared_ptr<ctgVariants> truth_vars = 
                clusterdata_ptr->ctg_superclusters[ctg]->ctg_variants[TRUTH][thi];
        int query_beg_idx = query_vars->clusters.size() ? query_vars->clusters[
            clusterdata_ptr->ctg_superclusters[ctg]->superclusters[QUERY][qhi][sc_idx]] : 0;
        int query_end_idx = query_vars->clusters.size() ? query_vars->clusters[
            clusterdata_ptr->ctg_superclusters[ctg]->superclusters[QUERY][qhi][sc_idx+1]] : 0;
        int truth_beg_idx = truth_vars->clusters.size() ? truth_vars->clusters[
            clusterdata_ptr->ctg_superclusters[ctg]->superclusters[TRUTH][thi][sc_idx]] : 0;
        int truth_end_idx = truth_vars->clusters.size() ? truth_vars->clusters[
            clusterdata_ptr->ctg_superclusters[ctg]->superclusters[TRUTH][thi][sc_idx+1]] : 0;

        // debug print
        if (print) printf("\nAlignment %s, aln_ptrs\n", aln_strs[i].data());
        if (print) printf("\nQUERY");
        if (print) print_ptrs(aln_ptrs, qi, query[i], truth[i]);
        if (print) printf("\nREF");
        if (print) print_ptrs(aln_ptrs, ri, ref, truth[i]);
        if (print) printf("Alignment %s, path_ptrs\n", aln_strs[i].data());
        if (print) printf("\nQUERY");
        if (print) print_ptrs(path_ptrs, qi, query[i], truth[i]);
        if (print) printf("\nREF");
        if (print) print_ptrs(path_ptrs, ri, ref, truth[i]);

        // init
        int ti_size = truth_ref_ptrs[i][PTRS].size();
        int hi = pr_query_ref_end[i];
        int qri_size = (hi % 2 == QUERY) ? query_ref_ptrs[i][PTRS].size() :
                ref_query_ptrs[i][PTRS].size();
        int prev_hi = hi;
        int qri = qri_size-1;
        int sync_ref_idx = ref_query_ptrs[i][PTRS].size();
        int prev_sync_ref_idx = sync_ref_idx;
        int prev_qri = qri;
        int ti = ti_size-1;
        int prev_ti = ti;
        int sync_truth_idx = ti_size;
        int prev_sync_truth_idx = ti_size;
        int new_ed = 0;
        int query_var_ptr = query_end_idx-1;
        int query_var_pos = query_vars->poss.size() && query_var_ptr >= 0 ? 
                query_vars->poss[query_var_ptr] - beg : 0;
        int prev_query_var_ptr = query_var_ptr;
        int truth_var_ptr = truth_end_idx-1;
        int truth_var_pos = truth_vars->poss.size() && truth_var_ptr >= 0 ? 
                truth_vars->poss[truth_var_ptr] - beg : 0;
        int prev_truth_var_ptr = truth_var_ptr;
        int path_idx = path[i].size()-1;

        if (print) printf("%s Path:\n", aln_strs[i].data());
        if (print) for (int p = path_idx; p >= 0; p--) {
            printf("(%s,%d,%d)\n", path[i][p].hi % 2 ? "REF" : "QRY", path[i][p].qri, path[i][p].ti);
        }

        if (print) printf("\n%s: %d\n", aln_strs[i].data(), beg);
        while (path_idx >= 0) {

            // set corresponding REF position
            int query_ref_pos = 0;
            if (hi == ri) { // REF, no variants
                query_ref_pos = qri;
            } else { // hi == qi, QUERY
                query_ref_pos = query_ref_ptrs[i][PTRS][qri];

                // in variant, ignore if didn't just finish REF FP or QRY PP/TP variant
                if (query_ref_ptrs[i][FLAGS][qri] & PTR_VARIANT) {
                    if(prev_hi == ri && !(ref_query_ptrs[i][FLAGS][prev_qri] & PTR_VAR_BEG)) {
                        query_ref_pos = query_var_pos+1; // skip next while loop
                    }
                    if (prev_hi == qi && (prev_qri == qri || 
                            !(query_ref_ptrs[i][FLAGS][prev_qri] & PTR_VAR_BEG))) { // not query var
                        query_ref_pos = query_var_pos+1; // skip next while loop
                    }
                }
            }
            // mark FPs and record variants within sync group
            while (query_ref_pos < query_var_pos && query_var_ptr >= query_beg_idx) { // just passed variant(s)
                if (prev_hi == ri) { // FP if when in variant was on REF
                    query_vars->errtypes[query_var_ptr] = ERRTYPE_FP;
                    query_vars->credit[query_var_ptr] = 0;
                    query_vars->callq[query_var_ptr] = 
                        query_vars->var_quals[query_var_ptr];
                    if (print) printf("QUERY REF='%s'\tALT='%s'\t%s\t%f\n",
                            query_vars->refs[query_var_ptr].data(),
                            query_vars->alts[query_var_ptr].data(), "FP", 0.0f);
                    query_var_ptr--;
                } else { // TP/PP
                    query_var_ptr--;
                }
                query_var_pos = (query_var_ptr < query_beg_idx) ? -1 :
                        query_vars->poss[query_var_ptr] - beg;
            }

            int truth_ref_pos = truth_ref_ptrs[i][PTRS][ti];
            if (truth_ref_ptrs[i][FLAGS][ti] & PTR_VARIANT) {
                if (prev_ti == ti || !(truth_ref_ptrs[i][FLAGS][prev_ti] & PTR_VAR_BEG))
                    truth_ref_pos = truth_var_pos + 1;
            }
            while (truth_ref_pos < truth_var_pos && truth_var_ptr >= truth_beg_idx) { // passed REF variant
                truth_var_ptr--;
                truth_var_pos = (truth_var_ptr < truth_beg_idx) ? -1 :
                    truth_vars->poss[truth_var_ptr] - beg;
            }

            // sync point: set TP/PP
            if (sync[i][path_idx]) {

                // calculate old edit distance
                int old_ed = 0;
                sync_ref_idx = (hi == ri) ? qri+1 : query_ref_ptrs[i][PTRS][qri]+1;
                sync_truth_idx = ti+1;
                std::vector< std::vector<int> > offs, ptrs;
                if (print) {
                    printf("  ref: %2d %s\n", int(ref.size()), ref.data());
                    printf("truth: %2d %s\n", int(truth[i].size()), truth[i].data());
                    printf("  ref[%2d:+%2d] =\n", sync_ref_idx, prev_sync_ref_idx-sync_ref_idx);
                    printf("%s\n", ref.substr(sync_ref_idx, 
                                prev_sync_ref_idx - sync_ref_idx).data());
                    printf("truth[%2d:+%2d] =\n", sync_truth_idx, prev_sync_truth_idx-sync_truth_idx);
                    printf("%s\n", truth[i].substr(sync_truth_idx, 
                                prev_sync_ref_idx - sync_ref_idx).data());
                }
                wf_ed(ref.substr(sync_ref_idx, 
                            prev_sync_ref_idx - sync_ref_idx), 
                        truth[i].substr(sync_truth_idx, 
                            prev_sync_truth_idx - sync_truth_idx), 
                        old_ed, offs, ptrs);

                if (old_ed == 0 && truth_var_ptr != prev_truth_var_ptr) 
                    ERROR("Old edit distance 0, variants exist.");

                // get min query var qual in sync section (for truth/query)
                float callq = g.max_qual;
                for (int query_var_idx = prev_query_var_ptr; 
                        query_var_idx > query_var_ptr; query_var_idx--) {
                    callq = std::min(callq, query_vars->var_quals[query_var_idx]);
                }

                // process QUERY variants
                for (int query_var_idx = prev_query_var_ptr; 
                        query_var_idx > query_var_ptr; query_var_idx--) {
                    float credit = 1 - float(new_ed)/old_ed;
                    // don't overwrite FPs
                    if (query_vars->errtypes[query_var_idx] == ERRTYPE_UN) {
                        if (new_ed == 0) { // TP
                            query_vars->errtypes[query_var_idx] = ERRTYPE_TP;
                            query_vars->credit[query_var_idx] = credit;
                            query_vars->callq[query_var_idx] = callq;
                            if (print) printf("QUERY REF='%s'\tALT='%s'\t%s\t%f\n",
                                    query_vars->refs[query_var_idx].data(),
                                    query_vars->alts[query_var_idx].data(), 
                                    "TP", credit);
                        } else if (new_ed == old_ed) { // FP
                            query_vars->errtypes[query_var_idx] = ERRTYPE_FP;
                            query_vars->credit[query_var_idx] = 0;
                            query_vars->callq[query_var_idx] = callq;
                            if (print) printf("QUERY REF='%s'\tALT='%s'\t%s\t%f\n",
                                    query_vars->refs[query_var_idx].data(),
                                    query_vars->alts[query_var_idx].data(), 
                                    "FP", 0.0f);
                        } else { // PP
                            query_vars->errtypes[query_var_idx] = ERRTYPE_PP;
                            query_vars->credit[query_var_idx] = credit;
                            query_vars->callq[query_var_idx] = callq;
                            if (print) printf("QUERY REF='%s'\tALT='%s'\t%s\t%f\n",
                                    query_vars->refs[query_var_idx].data(),
                                    query_vars->alts[query_var_idx].data(), 
                                    "PP", credit);
                        }
                    }
                }

                // process TRUTH variants
                for (int truth_var_idx = prev_truth_var_ptr; 
                        truth_var_idx > truth_var_ptr; truth_var_idx--) {
                    float credit = 1 - float(new_ed)/old_ed;
                    if (new_ed == 0) { // TP
                        truth_vars->errtypes[truth_var_idx] = ERRTYPE_TP;
                        truth_vars->credit[truth_var_idx] = credit;
                        truth_vars->callq[truth_var_idx] = callq;
                        if (print) printf("TRUTH REF='%s'\tALT='%s'\t%s\t%f\n",
                                truth_vars->refs[truth_var_idx].data(),
                                truth_vars->alts[truth_var_idx].data(), 
                                "TP", credit);
                    } else if (new_ed == old_ed) { // FP call, FN truth
                        truth_vars->errtypes[truth_var_idx] = ERRTYPE_FN;
                        truth_vars->credit[truth_var_idx] = credit;
                        truth_vars->callq[truth_var_idx] = g.max_qual;
                        if (print) printf("TRUTH REF='%s'\tALT='%s'\t%s\t%f\n",
                                truth_vars->refs[truth_var_idx].data(),
                                truth_vars->alts[truth_var_idx].data(), 
                                "FN", credit);
                    } else { // PP
                        truth_vars->errtypes[truth_var_idx] = ERRTYPE_PP;
                        truth_vars->credit[truth_var_idx] = credit;
                        truth_vars->callq[truth_var_idx] = callq;
                        if (print) printf("TRUTH REF='%s'\tALT='%s'\t%s\t%f\n",
                                truth_vars->refs[truth_var_idx].data(),
                                truth_vars->alts[truth_var_idx].data(), 
                                "PP", credit);
                    }
                }

                if (print) printf("SYNC @%s(%2d,%2d) = REF(%2d,%2d), ED %d->%d, QUERY %d-%d, TRUTH %d-%d\n",
                        (hi == ri) ? "REF" : "QRY", qri, ti,
                        (hi == ri) ? qri : query_ref_ptrs[i][PTRS][qri],
                        truth_ref_ptrs[i][PTRS][ti], old_ed, new_ed, 
                        prev_query_var_ptr, query_var_ptr,
                        prev_truth_var_ptr, truth_var_ptr
                );

                prev_query_var_ptr = query_var_ptr;
                prev_truth_var_ptr = truth_var_ptr;
                prev_sync_ref_idx = sync_ref_idx;
                prev_sync_truth_idx = sync_truth_idx;
                new_ed = 0;
            }

            // update pointers and edit distance
            if (print) printf("%s %s (%2d,%2d) %s\n", 
                    sync[i][path_idx] ? "*" : " ",
                    (hi == ri) ? "REF" : "QRY", qri, ti,
                    edits[i][path_idx] ? "X" : "=");

            // update path pointer
            path_idx--;
            if (path_idx < 0) break;

            // update location
            prev_qri = qri;
            qri = path[i][path_idx].qri;
            prev_ti = ti;
            ti = path[i][path_idx].ti;
            prev_hi = hi;
            hi = path[i][path_idx].hi;
            new_ed += edits[i][path_idx];

            // update next variant position
            query_var_pos = (query_var_ptr < query_beg_idx) ? -1 :
                query_vars->poss[query_var_ptr] - beg;
            truth_var_pos = (truth_var_ptr < truth_beg_idx) ? -1 :
                truth_vars->poss[truth_var_ptr] - beg;
        }
    }
}


/******************************************************************************/

void wf_ed(
        const std::string & query, const std::string & truth, int & s, 
        std::vector< std::vector<int> > & offs,
        std::vector< std::vector<int> > & ptrs
        ) {

    // alignment
    int query_len = query.size();
    int truth_len = truth.size();

    // early exit if either string is empty
    if (!query_len) { s = truth_len; return; }
    if (!truth_len) { s = query_len; return; }

    int mat_len = query_len + truth_len - 1;
    offs.push_back(std::vector<int>(mat_len, -2));
    offs[0][query_len-1] = -1;
    ptrs.push_back(std::vector<int>(mat_len, PTR_NONE));
    bool done = false;
    while (true) {

        // EXTEND WAVEFRONT
        for (int d = 0; d < mat_len; d++) {
            int off = offs[s][d];
            int diag = d + 1 - query_len;

            // don't allow starting from untouched cells
            if (off == -2) continue;

            // check that it's within matrix
            if (diag + off + 1 < 0) continue;
            if (off > query_len - 1) continue;
            if (diag + off > truth_len - 1) continue;

            // extend
            while (off < query_len - 1 && 
                   diag + off < truth_len - 1) {
                if (query[off+1] == truth[diag+off+1]) off++;
                else break;
            }
            offs[s][d] = off;

            // finish if done
            if (off == query_len - 1 && off + diag == truth_len - 1)
            { done = true; break; }

        }
        if (done) break;


        // NEXT WAVEFRONT
        // add wavefront, fill edge cells
        offs.push_back(std::vector<int>(mat_len, -2));
        ptrs.push_back(std::vector<int>(mat_len));
        // left cells
        if (query_len > 1 && s+1 == query_len-1) {
            offs[s+1][0] = s+1;
            ptrs[s+1][0] = PTR_INS;
        }
        // right cells
        if (truth_len > 1 && s+1 == mat_len-1) {
            offs[s+1][mat_len-1] = 0;
            ptrs[s+1][mat_len-1] = PTR_DEL;
        }

        // central cells
        for (int d = 1; d < mat_len-1; d++) {
            int offleft = offs[s][d-1];
            int offtop  = (offs[s][d+1] == -2) ? 
                -2 : offs[s][d+1]+1;
            int offdiag = (offs[s][d] == -2) ? 
                -2 : offs[s][d]+1;
            if (offdiag >= offtop && offdiag >= offleft) {
                offs[s+1][d] = offdiag;
                ptrs[s+1][d] = PTR_SUB;
            } else if (offleft >= offtop) {
                offs[s+1][d] = offleft;
                ptrs[s+1][d] = PTR_DEL;
            } else {
                offs[s+1][d] = offtop;
                ptrs[s+1][d] = PTR_INS;
            }
        }

        // single cell
        if (query_len == 1 && truth_len == 1) {
            offs[s+1][0] = offs[s][0]+1;
            ptrs[s+1][0] = PTR_SUB;
        }
        ++s;
    }
}

/******************************************************************************/

void wf_swg_align(
        const std::string & query, const std::string & truth, 
        std::vector< std::vector< std::vector<int> > > & ptrs,
        std::vector< std::vector< std::vector<int> > > & offs,
        int & s, int x, int o, int e, bool print
        ) {

    // init
    int query_len = query.size();
    int truth_len = truth.size();
    int mat_len = query_len + truth_len - 1;
    bool done = false;
    for (int m = 0; m < MATS; m++) {
        offs[m].push_back(std::vector<int>(mat_len, -2));
        ptrs[m].push_back(std::vector<int>(mat_len, PTR_NONE));
    }
    s = 0;
    offs[MAT_SUB][s][query_len-1] = -1; // main diag
    ptrs[MAT_SUB][s][query_len-1] = PTR_MAT;

    while (true) {

        // EXTEND WAVEFRONT (leave INS, DEL)
        for (int m = MAT_INS; m < MATS; m++) {
            for (int d = 0; d < mat_len; d++) {
                int off = offs[m][s][d];
                int diag = d + 1 - query_len;

                if (off >= 0 && off < query_len &&
                        diag+off >= 0 && diag+off < truth_len &&
                        offs[m][s][d] >= offs[MAT_SUB][s][d]) {
                    offs[MAT_SUB][s][d] = offs[m][s][d];
                    ptrs[MAT_SUB][s][d] |= (m == MAT_INS) ? PTR_INS : PTR_DEL;
                    if(print) printf("(S, %d, %d) swap\n", offs[MAT_SUB][s][d], 
                            offs[MAT_SUB][s][d]+d+1-query_len);

                }
            }
        }

        // EXTEND WAVEFRONT (diag, SUB only)
        for (int d = 0; d < mat_len; d++) {
            int off = offs[MAT_SUB][s][d];
            int diag = d + 1 - query_len;

            // extend
            while ( off != -2 && diag + off >= -1 && 
                    off < query_len - 1 && 
                    diag + off < truth_len - 1) {
                if (query[off+1] == truth[diag+off+1]) off++;
                else break;
            }
            if (off > offs[MAT_SUB][s][d])
                if(print) printf("(S, %d, %d) extend\n", off, off+diag);
            offs[MAT_SUB][s][d] = off;

            // finish if done
            if (off == query_len - 1 && off + diag == truth_len - 1)
            { done = true; break; }
        }
        if (done) break;

        for (int mi = 0; mi < MATS; mi++) {
            if(print) printf("\n%s matrix\n", type_strs[mi+1].data());
            if(print) printf("offs %d:", s);
            for (int di = 0; di < int(query.size() + truth.size()-1); di++) {
                if(print) printf("\t%d", offs[mi][s][di]);
            }
            if(print) printf("\n");
        }

        // NEXT WAVEFRONT
        s++;
        if(print) printf("\nscore = %d\n", s);
        for (int m = 0; m < MATS; m++) {
            offs[m].push_back(std::vector<int>(mat_len, -2));
            ptrs[m].push_back(std::vector<int>(mat_len, PTR_NONE));
        }

        for (int d = 0; d < mat_len; d++) {
            int diag = d + 1 - query_len;

            // sub (in SUB)
            if (s-x >= 0 && offs[MAT_SUB][s-x][d] != -2 && 
                            offs[MAT_SUB][s-x][d]+1 < query_len &&
                     diag + offs[MAT_SUB][s-x][d]+1 < truth_len &&
                            offs[MAT_SUB][s-x][d]+1 >= offs[MAT_SUB][s][d]) {
                offs[MAT_SUB][s][d] = offs[MAT_SUB][s-x][d] + 1;
                ptrs[MAT_SUB][s][d] |= PTR_SUB;
                if(print) printf("(S, %d, %d) sub\n", offs[MAT_SUB][s][d], 
                        offs[MAT_SUB][s][d]+diag);
            }

            // open gap (enter DEL)
            if (s-(o+e) >= 0 && d > 0 && 
                           offs[MAT_SUB][s-(o+e)][d-1] != -2 &&
                    diag + offs[MAT_SUB][s-(o+e)][d-1] < truth_len &&
                           offs[MAT_SUB][s-(o+e)][d-1] >= offs[MAT_DEL][s][d]) {
                offs[MAT_DEL][s][d] = offs[MAT_SUB][s-(o+e)][d-1];
                ptrs[MAT_DEL][s][d] |= PTR_SUB;
                if(print) printf("(D, %d, %d) open\n", offs[MAT_DEL][s][d], 
                        offs[MAT_DEL][s][d]+diag);
            }

            // open gap (enter INS)
            if (s-(o+e) >= 0 && d < mat_len-1 && 
                           offs[MAT_SUB][s-(o+e)][d+1] != -2 &&
                           offs[MAT_SUB][s-(o+e)][d+1]+1 < query_len &&
                    diag + offs[MAT_SUB][s-(o+e)][d+1]+1 < truth_len &&
                    diag + offs[MAT_SUB][s-(o+e)][d+1]+1 >= 0 &&
                           offs[MAT_SUB][s-(o+e)][d+1]+1 >= offs[MAT_INS][s][d]) {
                offs[MAT_INS][s][d] = offs[MAT_SUB][s-(o+e)][d+1]+1;
                ptrs[MAT_INS][s][d] |= PTR_SUB;
                if(print) printf("(I, %d, %d) open\n", offs[MAT_INS][s][d], 
                        offs[MAT_INS][s][d]+diag);
            }

            // extend gap (stay DEL)
            if (s-e >= 0 && d > 0 && 
                           offs[MAT_DEL][s-e][d-1] != -2 &&
                    diag + offs[MAT_DEL][s-e][d-1] < truth_len &&
                           offs[MAT_DEL][s-e][d-1] >= offs[MAT_DEL][s][d]) {
                offs[MAT_DEL][s][d] = offs[MAT_DEL][s-e][d-1];
                ptrs[MAT_DEL][s][d] |= PTR_DEL;
                if(print) printf("(D, %d, %d) extend\n", offs[MAT_DEL][s][d], 
                        offs[MAT_DEL][s][d]+diag);
            }

            // extend gap (stay INS)
            if (s-e >= 0 && d < mat_len-1 && 
                           offs[MAT_INS][s-e][d+1] != -2 &&
                           offs[MAT_INS][s-e][d+1]+1 < query_len &&
                    diag + offs[MAT_INS][s-e][d+1]+1 < truth_len &&
                    diag + offs[MAT_INS][s-e][d+1]+1 >= 0 &&
                           offs[MAT_INS][s-e][d+1]+1 >= offs[MAT_INS][s][d]) {
                offs[MAT_INS][s][d] = offs[MAT_INS][s-e][d+1]+1;
                ptrs[MAT_INS][s][d] |= PTR_INS;
                if(print) printf("(I, %d, %d) extend\n", offs[MAT_INS][s][d], 
                        offs[MAT_INS][s][d]+diag);
            }
        }
    }
}

/******************************************************************************/


editData alignment_wrapper(std::shared_ptr<superclusterData> clusterdata_ptr) {
    if (g.verbosity >= 1) INFO(" ");
    if (g.verbosity >= 1) INFO("Calculating precision/recall and distance metrics");

    // +2 since it's inclusive, but then also needs to include one quality higher
    // which doesn't contain any variants (to get draft reference edit dist)
    std::vector<int> all_qual_dists(g.max_qual+2, 0);
    editData edits;
    for (std::string ctg : clusterdata_ptr->contigs) {
        std::vector<int> ctg_qual_dists(g.max_qual+2,0);
        if (clusterdata_ptr->ctg_superclusters[ctg]->n && g.verbosity >= 1)
            INFO("  Contig '%s'", ctg.data())

        // set superclusters pointer
        std::shared_ptr<ctgSuperclusters> sc = clusterdata_ptr->ctg_superclusters[ctg];

        // iterate over superclusters
        for(int sc_idx = 0; sc_idx < clusterdata_ptr->ctg_superclusters[ctg]->n; sc_idx++) {

            /////////////////////////////////////////////////////////////////////
            // PRECISION-RECALL: allow skipping called variants                  
            /////////////////////////////////////////////////////////////////////
            
            // set pointers between each hap (query1/2, truth1/2) and reference
            std::string query1 = "", ref_q1 = ""; 
            std::vector< std::vector<int> > query1_ref_ptrs, ref_query1_ptrs;
            generate_ptrs_strs(
                    query1, ref_q1, query1_ref_ptrs, ref_query1_ptrs, 
                    sc->ctg_variants[QUERY][HAP1], 
                    sc->superclusters[QUERY][HAP1][sc_idx],
                    sc->superclusters[QUERY][HAP1][sc_idx+1],
                    sc->begs[sc_idx], sc->ends[sc_idx], clusterdata_ptr->ref, ctg
            );
            std::string query2 = "", ref_q2 = ""; 
            std::vector< std::vector<int> > query2_ref_ptrs, ref_query2_ptrs;
            generate_ptrs_strs(
                    query2, ref_q2, query2_ref_ptrs, ref_query2_ptrs, 
                    sc->ctg_variants[QUERY][HAP2],
                    sc->superclusters[QUERY][HAP2][sc_idx],
                    sc->superclusters[QUERY][HAP2][sc_idx+1],
                    sc->begs[sc_idx], sc->ends[sc_idx], clusterdata_ptr->ref, ctg
            );
            std::string truth1 = "", ref_t1 = ""; 
            std::vector< std::vector<int> > truth1_ref_ptrs, ref_truth1_ptrs;
            generate_ptrs_strs(
                    truth1, ref_t1, truth1_ref_ptrs, ref_truth1_ptrs, 
                    sc->ctg_variants[TRUTH][HAP1],
                    sc->superclusters[TRUTH][HAP1][sc_idx],
                    sc->superclusters[TRUTH][HAP1][sc_idx+1],
                    sc->begs[sc_idx], sc->ends[sc_idx], clusterdata_ptr->ref, ctg
            );
            std::string truth2 = "", ref_t2 = ""; 
            std::vector< std::vector<int> > truth2_ref_ptrs, ref_truth2_ptrs;
            generate_ptrs_strs(
                    truth2, ref_t2, truth2_ref_ptrs, ref_truth2_ptrs, 
                    sc->ctg_variants[TRUTH][HAP2],
                    sc->superclusters[TRUTH][HAP2][sc_idx],
                    sc->superclusters[TRUTH][HAP2][sc_idx+1],
                    sc->begs[sc_idx], sc->ends[sc_idx], clusterdata_ptr->ref, ctg
            );

            if (false) {
                printf("\n%s:%d\n", ctg.data(), sc->begs[sc_idx]);
                printf("REF:       %s\n", ref_q1.data());
                printf("QUERY1:    %s\n", query1.data());
                printf("QUERY2:    %s\n", query2.data());
                printf("TRUTH1:    %s\n", truth1.data());
                printf("TRUTH2:    %s\n", truth2.data());
                printf("Q1->REF PTRS: "); print_ref_ptrs(query1_ref_ptrs);
                printf("REF->Q1 PTRS: "); print_ref_ptrs(ref_query1_ptrs);
                printf("Q2->REF PTRS: "); print_ref_ptrs(query2_ref_ptrs);
                printf("REF->Q2 PTRS: "); print_ref_ptrs(ref_query2_ptrs);
                printf("T1->REF PTRS: "); print_ref_ptrs(truth1_ref_ptrs);
                printf("REF->T1 PTRS: "); print_ref_ptrs(ref_truth1_ptrs);
                printf("T2->REF PTRS: "); print_ref_ptrs(truth2_ref_ptrs);
                printf("REF->T2 PTRS: "); print_ref_ptrs(ref_truth2_ptrs);
            }

            // calculate four forward-pass alignment edit dists
            // query1-truth2, query1-truth1, query2-truth1, query2-truth2
g.timers[TIME_PR].start();
            std::vector<int> aln_score(4), aln_query_ref_end(4);
            std::unordered_map<idx1, int> aln_ptrs;
            std::vector< std::unordered_map<idx1, idx1> > swap_pred_map(4);
            calc_prec_recall_aln(
                    query1, query2, truth1, truth2, ref_q1,
                    query1_ref_ptrs, ref_query1_ptrs, 
                    query2_ref_ptrs, ref_query2_ptrs,
                    truth1_ref_ptrs, truth2_ref_ptrs,
                    aln_score, aln_ptrs, swap_pred_map,
                    aln_query_ref_end, false
            );

            // store optimal phasing for each supercluster
            // ORIG: query1-truth1 and query2-truth2
            // SWAP: query1-truth2 and query2-truth1
            int phase = store_phase(clusterdata_ptr, ctg, aln_score);

            // calculate paths from alignment
            std::vector< std::vector<idx1> > path;
            std::vector< std::vector<bool> > sync, edit;
            std::unordered_map<idx1, int> path_ptrs, path_scores;
            calc_prec_recall_path(
                    ref_q1, query1, query2, truth1, truth2,
                    path, sync, edit,
                    aln_ptrs, path_ptrs, path_scores,
                    query1_ref_ptrs, ref_query1_ptrs, 
                    query2_ref_ptrs, ref_query2_ptrs, 
                    truth1_ref_ptrs, truth2_ref_ptrs,
                    swap_pred_map,
                    aln_query_ref_end, false);

            // calculate precision/recall from paths
            calc_prec_recall(
                    clusterdata_ptr, sc_idx, 
                    ctg, 
                    ref_q1,
                    query1, query2, truth1, truth2,
                    path, sync, edit,
                    aln_ptrs, path_ptrs, 
                    query1_ref_ptrs, ref_query1_ptrs, 
                    query2_ref_ptrs, ref_query2_ptrs,
                    truth1_ref_ptrs, truth2_ref_ptrs,
                    aln_query_ref_end, phase, 
                    false
            );
g.timers[TIME_PR].stop();


            /////////////////////////////////////////////////////////////////////
            // SMITH-WATERMAN DISTANCE: don't allow skipping called variants     
            /////////////////////////////////////////////////////////////////////
            
            // keep or swap truth haps based on previously decided phasing
            std::vector<std::string> truth(2);
            if (phase == PHASE_SWAP) {
                truth[HAP1] = truth2; truth[HAP2] = truth1;
            } else {
                truth[HAP1] = truth1; truth[HAP2] = truth2;
            }

            // phasing is known, add scores for each hap
            std::vector<int> qual_dists(g.max_qual+2, 0);
            for (int hap = 0; hap < HAPS; hap++) {

                // supercluster start/end indices
                int beg_idx = sc->ctg_variants[QUERY][hap]->clusters.size() ?
                    sc->ctg_variants[QUERY][hap]->clusters[
                        sc->superclusters[QUERY][hap][sc_idx]] : 0;
                int end_idx = sc->ctg_variants[QUERY][hap]->clusters.size() ?
                    sc->ctg_variants[QUERY][hap]->clusters[
                        sc->superclusters[QUERY][hap][sc_idx+1]] : 0;

                // calculate quality thresholds 
                // (where string would change when including/excluding variants)
                std::set<int> quals = {};
                for (int var_idx = beg_idx; var_idx < end_idx; var_idx++) {
                    quals.insert(sc->ctg_variants[QUERY][hap]->var_quals[var_idx]+1);
                }
                quals.insert(g.max_qual+2);

                // sweep through quality thresholds
                int prev_qual = 0;
                for (int qual : quals) {
                    if(false) printf("qual: %d-%d\n", prev_qual, qual-1);

                    // generate query string (only applying variants with Q>=qual)
                    std::string query = generate_str(
                            clusterdata_ptr->ref, 
                            sc->ctg_variants[QUERY][hap], 
                            ctg, beg_idx, end_idx,
                            sc->begs[sc_idx], sc->ends[sc_idx], 
                            prev_qual);

                    // align strings, backtrack, calculate distance
g.timers[TIME_SW].start();
                    std::vector< std::vector< std::vector<int> > > ptrs(MATS);
                    std::vector< std::vector< std::vector<int> > > offs(MATS);
                    int s = 0;
                    std::reverse(query.begin(), query.end());
                    std::reverse(truth[hap].begin(), truth[hap].end());
                    wf_swg_align(query, truth[hap], ptrs, offs,
                            s, g.eval_sub, g.eval_open, g.eval_extend, false);
g.timers[TIME_SW].stop();
                    std::vector<int> cigar = wf_swg_backtrack(query, truth[hap], 
                            ptrs, offs, s, g.eval_sub, g.eval_open, g.eval_extend, false);
                    std::reverse(query.begin(), query.end());
                    std::reverse(truth[hap].begin(), truth[hap].end());
                    std::reverse(cigar.begin(), cigar.end());
                    int dist = count_dist(cigar);

                    // add distance for range of corresponding quals
                    for (int q = prev_qual; q < qual; q++) {
                        all_qual_dists[q] += dist;
                        ctg_qual_dists[q] += dist;
                        qual_dists[q] += dist;
                        edits.add_edits(ctg, sc->begs[sc_idx], hap, cigar, sc_idx, q);
                    }
                    prev_qual = qual;
                }
            }


            /////////////////////////////////////////////////////////////////////
            // DEBUG PRINTING                                                    
            /////////////////////////////////////////////////////////////////////
            if (false) {
                // print cluster info
                printf("\n\nSupercluster: %d\n", sc_idx);
                for (int i = 0; i < CALLSETS*HAPS; i++) {
                    int callset = i >> 1;
                    int hap = i % 2;
                    int cluster_beg = sc->superclusters[callset][hap][sc_idx];
                    int cluster_end = sc->superclusters[callset][hap][sc_idx+1];
                    printf("%s%d: %d clusters (%d-%d)\n", 
                        callset_strs[callset].data(), hap+1,
                        cluster_end-cluster_beg,
                        cluster_beg, cluster_end);

                    for (int j = cluster_beg; j < cluster_end; j++) {
                        auto vars = sc->ctg_variants[callset][hap];
                        int variant_beg = vars->clusters[j];
                        int variant_end = vars->clusters[j+1];
                        printf("\tCluster %d: %d variants (%d-%d)\n", j, 
                            variant_end-variant_beg, variant_beg, variant_end);
                        for (int k = variant_beg; k < variant_end; k++) {
                            printf("\t\t%s %d\t%s\t%s\tQ=%f\n", ctg.data(), vars->poss[k], 
                            vars->refs[k].size() ?  vars->refs[k].data() : "_", 
                            vars->alts[k].size() ?  vars->alts[k].data() : "_",
                            vars->var_quals[k]);
                        }
                    }
                }
                printf("Edit Distance: %d\n", 
                    *std::min_element(qual_dists.begin(), qual_dists.end()));

            } // debug print

        } // each cluster
        if (clusterdata_ptr->ctg_superclusters[ctg]->n && g.verbosity >= 1)
            INFO("    %d edits", 
                *std::min_element(ctg_qual_dists.begin(), ctg_qual_dists.end()));
    } // each contig
    if (g.verbosity >= 1) INFO(" ");
    if (g.verbosity >= 1) INFO("  Total edit distance: %d", 
                *std::min_element(all_qual_dists.begin(), all_qual_dists.end()));
    return edits;
}


/******************************************************************************/


/* Calculate the combined Smith-Waterman score for all variants within 
 * one or several adjacent clusters.
 */
int calc_vcf_swg_score(std::shared_ptr<ctgVariants> vars, 
        int clust_beg_idx, int clust_end_idx, int sub, int open, int extend) {
    int score = 0;
    for (int var_idx = vars->clusters[clust_beg_idx]; 
            var_idx < vars->clusters[clust_end_idx]; var_idx++) {
        switch (vars->types[var_idx]) {
            case TYPE_SUB:
                score += sub;
                break;
            case TYPE_INS:
                score += open;
                score += extend * vars->alts[var_idx].size();
                break;
            case TYPE_DEL:
                score += open;
                score += extend * vars->refs[var_idx].size();
                break;
            default:
                ERROR("Unexpected variant type in calc_vcf_swg_score()");
        }
    }
    return score;
}


int calc_cig_swg_score(const std::vector<int> & cigar, 
        int sub, int open, int extend) {
    int score = 0;
    int prev = PTR_MAT;
    for (int i = 0; i < int(cigar.size()); i++) {
        switch (cigar[i]) {
            case PTR_MAT:
                i++; // increment by 2
                prev = PTR_MAT;
                break;
            case PTR_SUB:
                score += sub;
                prev = PTR_SUB;
                i++; // increment by 2
                break;
            case PTR_INS:
                score += (prev == PTR_INS ? extend : open + extend);
                prev = PTR_INS;
                break;
            case PTR_DEL:
                score += (prev == PTR_DEL ? extend : open + extend);
                prev = PTR_DEL;
                break;
            default:
                ERROR("Unexpected variant type in calc_cig_swg_score()");
        }
    }
    return score;
}


/******************************************************************************/


/* Perform Djikstra Smith-Waterman alignment of two strings, returning 
 * the farthest-reaching reference index of lesser or equal score to that provided.
 * Strings are generated by applying variants to draft ref, and skipping 
 * variants is not allowed. Neither is a diagonal ref transition.
 */
int swg_max_reach(const std::string & query, const std::string & ref, 
        const std::vector< std::vector<int> > & query_ref_ptrs,
        const std::vector< std::vector<int> > & ref_query_ptrs,
        int sub, int open, int extend, int score, bool print,
        bool reverse /*= false*/, int ref_section /*= -1*/) {
    
    int ref_len = ref.size();
    int query_len = query.size();
    int sub_wave = 0, open_wave = 0, extend_wave = 0, end_wave = 0;
        
    // set first wavefront
    std::queue<idx2> queue;
    queue.push({MAT_SUB, 0, 0});
    std::unordered_set<idx2> done;
    done.insert({MAT_SUB, 0, 0});
    std::vector< std::unordered_set<idx2> > waves(score+1);

    // find furthest-reaching alignment with lesser or equal score
    for (int s = 0; s <= score; s++) {
        if (print) printf("s = %d\n", s);

        // EXTEND WAVEFRONT (stay at same score)
        while (!queue.empty()) {
            idx2 x = queue.front(); queue.pop();
            if (print) printf("  x = (%c, %d, %d)\n", 
                    std::string("SID")[x.mi], x.qi, x.ri);
            waves[s].insert(x);

            // allow non-main-diagonal match
            bool main_diag = x.ri >= ref_section &&
                    x.qi+1 < query_len && x.ri+1 < ref_len && // repeat bounds check
                    !(query_ref_ptrs[FLAGS][x.qi] & PTR_VARIANT) && // neither prev/curr are variants
                    !(ref_query_ptrs[FLAGS][x.ri] & PTR_VARIANT) &&
                    !(query_ref_ptrs[FLAGS][x.qi+1] & PTR_VARIANT) &&
                    !(ref_query_ptrs[FLAGS][x.ri+1] & PTR_VARIANT) &&
                    query_ref_ptrs[PTRS][x.qi] == x.ri && // point to one another (main diag)
                    ref_query_ptrs[PTRS][x.ri] == x.qi &&
                    query_ref_ptrs[PTRS][x.qi+1] == x.ri+1 &&
                    ref_query_ptrs[PTRS][x.ri+1] == x.qi+1;
            if (x.qi+1 < query_len && x.ri+1 < ref_len && 
                    !main_diag &&
                    query[x.qi+1] == ref[x.ri+1] && 
                    !contains(done, {x.mi, x.qi+1, x.ri+1})) {
                idx2 next(x.mi, x.qi+1, x.ri+1);
                queue.push(next);
                done.insert(next);
            }

            // allow exiting D/I state (no penalty for forward only)
            if (!reverse && (x.mi == MAT_INS || x.mi == MAT_DEL) && 
                    !contains(done, {MAT_SUB, x.qi, x.ri})) {
                queue.push({MAT_SUB, x.qi, x.ri});
                done.insert({MAT_SUB, x.qi, x.ri});
            }
        }

        // exit if we're done aligning
        if (s == score || contains(done, {MAT_SUB, query_len-1, ref_len-1})) 
            break;

        // NEXT WAVEFRONT (increase score by one)
        
        // SUB transition (only from MAT_SUB)
        sub_wave = s + 1 - sub;
        if (sub_wave >= 0 && sub_wave <= score) {
            for (idx2 x : waves[sub_wave]) {
                if (x.mi == MAT_SUB && x.qi+1 < query_len && x.ri+1 < ref_len) {
                    idx2 next({x.mi, x.qi+1, x.ri+1});
                    if (!contains(done, next)) {
                        queue.push(next);
                        done.insert(next);
                    }
                }
            }
        }

        // INS/DEL opening (only from MAT_SUB)
        if (reverse) // g.open penalty when leaving INDEL
            open_wave = s + 1 - extend;
        else // g.open penalty when entering INDEL
            open_wave = s + 1 - (open+extend);
        if (open_wave >= 0 && open_wave <= score) {
            for(idx2 x : waves[open_wave]) {

                // INS opening
                if (x.mi == MAT_SUB && x.qi+1 < query_len) {
                    idx2 next({MAT_INS, x.qi+1, x.ri});
                    if (!contains(done, next)) {
                        queue.push(next);
                        done.insert(next);
                    }
                }

                // DEL opening
                if (x.mi == MAT_SUB && x.ri+1 < ref_len) {
                    idx2 next({MAT_DEL, x.qi, x.ri+1});
                    if (!contains(done, next)) {
                        queue.push(next);
                        done.insert(next);
                    }
                }
            }
        }

        // reverse alignment, penalize leaving INDEL
        end_wave = s + 1 - open;
        if (reverse && end_wave >= 0 && end_wave <= score) {
            for (idx2 x : waves[end_wave]) {
                if ((x.mi == MAT_INS || x.mi == MAT_DEL) && 
                        !contains(done, {MAT_SUB, x.qi, x.ri})) {
                    queue.push({MAT_SUB, x.qi, x.ri});
                    done.insert({MAT_SUB, x.qi, x.ri});
                }
            }
        }

        // INS/DEL extension
        extend_wave = s + 1 - extend;
        if (extend_wave >= 0 && extend_wave <= score) {
            for(idx2 x : waves[extend_wave]) {

                // INS extending (only from MAT_INS)
                if (x.mi == MAT_INS && x.qi+1 < query_len) {
                    idx2 next({x.mi, x.qi+1, x.ri});
                    if (!contains(done, next)) {
                        queue.push(next);
                        done.insert(next);
                    }
                }

                // DEL extending (only from MAT_DEL)
                if (x.mi == MAT_DEL && x.ri+1 < ref_len) {
                    idx2 next({x.mi, x.qi, x.ri+1});
                    if (!contains(done, next)) {
                        queue.push(next);
                        done.insert(next);
                    }
                }
            }
        }
    } // end reach

    // farthest reach isn't necessarily at max score
    int min_wave = std::max(0, std::min(std::min(
                    open_wave, extend_wave), sub_wave));
    // get max reach
    int max_reach = 0;
    for (int w = min_wave; w <= score; w++) {
        for (idx2 x : waves[w]) {
            max_reach = std::max(max_reach, x.ri);
        }
    }
    return max_reach;
}


/******************************************************************************/


std::unique_ptr<variantData> wf_swg_realign(
        std::unique_ptr<variantData> & vcf, 
        std::shared_ptr<fastaData> ref_fasta, 
        int sub, int open, int extend, int callset, bool print /* = false */) {
    if (g.verbosity >= 1) INFO(" ");
    if (g.verbosity >= 1) INFO("Realigning %s VCF '%s'", 
            callset_strs[callset].data(), vcf->filename.data());

    // copy vcf header data over to results vcf
    std::unique_ptr<variantData> results(new variantData());
    results->set_header(vcf);

    // iterate over each contig haplotype
    for (int hap = 0; hap < 2; hap++) {
        for (auto itr = vcf->ctg_variants[hap].begin(); 
                itr != vcf->ctg_variants[hap].end(); itr++) {
            std::string ctg = itr->first;
            std::shared_ptr<ctgVariants> vars = itr->second;
            if (vars->poss.size() == 0) continue;
            if (g.verbosity >= 1) INFO("  Haplotype %d Contig %s", hap+1, ctg.data());

            // realign each cluster of variants
            for (int cluster = 0; cluster < int(vars->clusters.size()-1); cluster++) {
                int beg_idx = vars->clusters[cluster];
                int end_idx = vars->clusters[cluster+1];
                int beg = vars->poss[beg_idx]-1;
                int end = vars->poss[end_idx-1] + vars->rlens[end_idx-1]+1;

                // variant qual is minimum in cluster
                float qual = g.max_qual;
                for (int i = beg_idx; i < end_idx; i++) {
                    qual = std::min(qual, vars->var_quals[i]);
                }

                // generate strings
                std::string query = 
                    generate_str(ref_fasta, vars, ctg, beg_idx, end_idx, beg, end);
                std::string ref = ref_fasta->fasta.at(ctg).substr(beg, end-beg);
                
                // perform alignment
                if (print) printf("REF:   %s\n", ref.data());
                if (print) printf("QUERY: %s\n", query.data());
                std::vector< std::vector< std::vector<int> > > ptrs(MATS);
                std::vector< std::vector< std::vector<int> > > offs(MATS);
                int s = 0;
                std::reverse(query.begin(), query.end()); // for left-aligned INDELs
                std::reverse(ref.begin(), ref.end());
                wf_swg_align(query, ref, ptrs, offs, s, sub, open, extend, false);
                
                // backtrack
                std::vector<int> cigar = wf_swg_backtrack(query, ref, ptrs, offs, 
                        s, sub, open, extend, print);
                std::reverse(query.begin(), query.end());
                std::reverse(ref.begin(), ref.end());
                std::reverse(cigar.begin(), cigar.end());
                if (print) print_cigar(cigar);

                // compare distances
                if (print) {
                    int new_score = calc_cig_swg_score(cigar, sub, open, extend);
                    int old_score = calc_vcf_swg_score(vars, cluster, cluster+1, 
                            sub, open, extend);
                    if (new_score < old_score) {
                        printf("\n\tCluster %d: %d variants (%d-%d)\n", 
                            cluster, end_idx-beg_idx, beg_idx, end_idx);
                        for (int i = beg_idx; i < end_idx; i++) {
                            printf("\t\t%s %d\t%s\t%s\tQ=%f\n", 
                                ctg.data(), vars->poss[i], 
                                vars->refs[i].size() ?  vars->refs[i].data() : "_", 
                                vars->alts[i].size() ?  vars->alts[i].data() : "_",
                                vars->var_quals[i]);
                        }
                        printf("Old score: %d\n", old_score);
                        printf("New score: %d\n", new_score);
                    }
                }
                
                // save resulting variants
                results->add_variants(cigar, hap, beg, ctg, query, ref, qual);

            } // cluster
        } // contig
    } // hap

    return results;
}


/******************************************************************************/


int count_dist(const std::vector<int> & cigar) {
    int cigar_ptr = 0;
    int dist = 0;
    while (cigar_ptr < int(cigar.size())) {
        switch (cigar[cigar_ptr]) {
        case PTR_MAT:
            cigar_ptr += 2;
            break;
        case PTR_SUB:
            dist++;
            cigar_ptr += 2;
            break;
        case PTR_INS:
        case PTR_DEL:
            dist++;
            cigar_ptr++;
            break;
        default:
            ERROR("Unexpected CIGAR in count_dist()");
        }
    }
    return dist;
}


/******************************************************************************/


/* After alignment, backtrack using pointers and save cigar. */
std::vector<int> wf_swg_backtrack(
        const std::string & query,
        const std::string & ref,
        const std::vector< std::vector< std::vector<int> > > & ptrs, 
        const std::vector< std::vector< std::vector<int> > > & offs, 
        int s, int x, int o, int e,
        bool print) {

    /* if (print) print_wfa_ptrs(query, ref, s, ptrs, offs); */
    /* for (int mi = 0; mi < MATS; mi++) { */
    /*     if(print) printf("\n%s matrix\n", type_strs[mi+1].data()); */
    /*     for (int si = 0; si <= s; si++) { */
    /*         if(print) printf("%d:\n", si); */
    /*         for (int di = 0; di < int(query.size() + ref.size()-1); di++) { */
    /*             if(print) printf("\t%d", offs[mi][si][di]); */
    /*         } */
    /*         if(print) printf("\n"); */
    /*     } */
    /* } */

    // initialize backtrack position and cigar string
    std::vector<int> cigar(query.size() + ref.size(), 0);
    int cigar_ptr = cigar.size()-1;
    idx2 pos(MAT_SUB, query.size()-1, ref.size()-1);

    if (print) printf("\nWF SWG Backtrack:\n");
    while (pos.qi >= 0 || pos.ri >= 0) {

        // debug print
        if (s < 0) ERROR("Negative wf_swg_backtrack() score at (%c, %d, %d)",
                std::string("SID")[pos.mi], pos.qi, pos.ri);

        // init
        int diag = pos.ri - pos.qi;
        int d = query.size()-1 + diag;
        
        if (pos.mi == MAT_SUB) {
            if (ptrs[MAT_SUB][s][d] & PTR_INS) { // end INS freely
                int prev_off = offs[MAT_INS][s][d];
                while (pos.qi > prev_off) { // match
                    if (print) printf("(%c, %d, %d) match ins\n",
                            std::string("SID")[pos.mi], pos.qi, pos.ri);
                    pos.qi--; cigar[cigar_ptr--] = PTR_MAT;
                    pos.ri--; cigar[cigar_ptr--] = PTR_MAT;
                    if (pos.qi < 0 || pos.ri < 0) ERROR("During wf_swg_backtrack() at (%c, %d, %d)",
                            std::string("SID")[pos.mi], pos.qi, pos.ri);
                }
                pos.mi = MAT_INS;
            } else if (ptrs[MAT_SUB][s][d] & PTR_DEL) { // end DEL freely
                int prev_off = offs[MAT_DEL][s][d];
                while (pos.qi > prev_off) { // match
                    if (print) printf("(%c, %d, %d) match del\n",
                            std::string("SID")[pos.mi], pos.qi, pos.ri);
                    pos.qi--; cigar[cigar_ptr--] = PTR_MAT;
                    pos.ri--; cigar[cigar_ptr--] = PTR_MAT;
                    if (pos.qi < 0 || pos.ri < 0) ERROR("During wf_swg_backtrack() at (%c, %d, %d)",
                            std::string("SID")[pos.mi], pos.qi, pos.ri);
                }
                pos.mi = MAT_DEL;
            } else if (ptrs[MAT_SUB][s][d] & PTR_SUB) { // sub
                if (s-x < 0) ERROR("Unexpected PTR_SUB");
                int prev_off = offs[MAT_SUB][s-x][d];
                while (pos.qi > prev_off+1) { // match
                    if (print) printf("(%c, %d, %d) match sub\n",
                            std::string("SID")[pos.mi], pos.qi, pos.ri);
                    pos.qi--; cigar[cigar_ptr--] = PTR_MAT;
                    pos.ri--; cigar[cigar_ptr--] = PTR_MAT;
                    if (pos.qi < 0 || pos.ri < 0) ERROR("During wf_swg_backtrack() at (%c, %d, %d)",
                            std::string("SID")[pos.mi], pos.qi, pos.ri);
                }
                if (print) printf("(%c, %d, %d) sub\n",
                        std::string("SID")[pos.mi], pos.qi, pos.ri);
                pos.qi--; cigar[cigar_ptr--] = PTR_SUB;
                pos.ri--; cigar[cigar_ptr--] = PTR_SUB;
                s -= x;
            } else if (ptrs[MAT_SUB][s][d] & PTR_MAT) { // only matches remain
                while (pos.qi >= 0 && pos.ri >= 0) { // match
                    if (print) printf("(%c, %d, %d) match mat\n",
                            std::string("SID")[pos.mi], pos.qi, pos.ri);
                    pos.qi--; cigar[cigar_ptr--] = PTR_MAT;
                    pos.ri--; cigar[cigar_ptr--] = PTR_MAT;
                }
                if (pos.qi >= 0 || pos.ri >= 0)
                    ERROR("PTR_MAT expected all matches, now at (%c, %d, %d)",
                        std::string("SID")[pos.mi], pos.qi, pos.ri);
            } else {
                ERROR("Unexpected pointer '%d' in wf_swg_backtrack() at (%c, %d, %d)",
                    ptrs[MAT_SUB][s][d], std::string("SID")[pos.mi], pos.qi, pos.ri);
            }

        } else if (pos.mi == MAT_INS) {
            if (ptrs[MAT_INS][s][d] & PTR_INS) { // extend INS
                if (print) printf("(%c, %d, %d) ext ins\n",
                        std::string("SID")[pos.mi], pos.qi, pos.ri);
                pos.qi--; cigar[cigar_ptr--] = PTR_INS;
                s -= e;
            } else if (ptrs[MAT_INS][s][d] & PTR_SUB) { // start INS
                if (print) printf("(%c, %d, %d) start ins\n",
                        std::string("SID")[pos.mi], pos.qi, pos.ri);
                pos.qi--; cigar[cigar_ptr--] = PTR_INS;
                pos.mi = MAT_SUB;
                s -= o + e;
            } else {
                ERROR("Unexpected pointer '%d' in wf_swg_backtrack() at (%c, %d, %d)",
                    ptrs[MAT_INS][s][d], std::string("SID")[pos.mi], pos.qi, pos.ri);
            }

        } else if (pos.mi == MAT_DEL) {
            if (ptrs[MAT_DEL][s][d] & PTR_DEL) { // extend DEL
                if (print) printf("(%c, %d, %d) ext del\n",
                        std::string("SID")[pos.mi], pos.qi, pos.ri);
                pos.ri--; cigar[cigar_ptr--] = PTR_DEL;
                s -= e;
            } else if (ptrs[MAT_DEL][s][d] & PTR_SUB) { // start DEL
                if (print) printf("(%c, %d, %d) start del\n",
                        std::string("SID")[pos.mi], pos.qi, pos.ri);
                pos.ri--; cigar[cigar_ptr--] = PTR_DEL;
                pos.mi = MAT_SUB;
                s -= o + e;
            } else {
                ERROR("Unexpected pointer '%d' in wf_swg_backtrack() at (%c, %d, %d)",
                    ptrs[MAT_DEL][s][d], std::string("SID")[pos.mi], pos.qi, pos.ri);
            }

        } else {
            ERROR("Unexpected MAT type in wf_swg_backtrack()");
        }
        if (!(pos.qi == -1 && pos.ri == -1) && (pos.qi < 0 || pos.ri < 0)) 
            ERROR("During wf_swg_backtrack() at (%c, %d, %d)",
                std::string("SID")[pos.mi], pos.qi, pos.ri);
    }
    return cigar;
}
