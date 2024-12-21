#include <string>
#include <unordered_map>
#include <vector>
#include <cmath>

#include "htslib/vcf.h"

#include "variant.h"
#include "print.h"
#include "dist.h"

/******************************************************************************/

/* Merge unclustered variants with another set of unclustered variants. */
void variantData::merge(std::shared_ptr<variantData> other_variants) {
    for (const std::string & ctg : other_variants->contigs) {

        // if another contig is completely new, just store the whole thing
        if (std::find(this->contigs.begin(), this->contigs.end(), ctg) == this->contigs.end()) {
            for (int hi = 0; hi < HAPS; hi++) {
                this->variants[hi][ctg] = other_variants->variants[hi][ctg];
            }
            continue;
        }

        // otherwise, contig exists for both, add variants in correct order
        for (int hi = 0; hi < HAPS; hi++) {
            int vi1 = 0;
            int vi2 = 0;
            std::shared_ptr<ctgVariants> merged_vars(new ctgVariants(ctg));
            std::shared_ptr<ctgVariants> v1 = this->variants[hi][ctg];
            std::shared_ptr<ctgVariants> v2 = other_variants->variants[hi][ctg];
            while (vi1 < v1->n && vi2 < v2->n) {
                if (v1->poss[vi1] < v2->poss[vi2]) { // v1 occurs first
                    merged_vars->add_var(v1, vi1++);
                } else if (v1->poss[vi1] > v2->poss[vi2]) { // v2 occurs first
                    merged_vars->add_var(v2, vi2++);
                } else { // tie, order by variant type
                    // TODO: what if variants are the same ref/alt?
                    if (v1->types[vi1] == TYPE_INS) {
                        merged_vars->add_var(v1, vi1++);
                    } else if (v2->types[vi2] == TYPE_INS) {
                        merged_vars->add_var(v2, vi2++);
                    } else {
                        merged_vars->add_var(v1, vi1++);
                    }
                }
            }
            // add all remaining vars
            while (vi1 < v1->n) { merged_vars->add_var(v1, vi1++); }
            while (vi2 < v2->n) { merged_vars->add_var(v2, vi2++); }
            this->variants[hi][ctg] = merged_vars;
        }
    }
}

/******************************************************************************/

ctgVariants::ctgVariants(const std::string & ctg) { 
    this->ctg = ctg;
    this->n = 0; 
    for (int i = 0; i < PHASES; i++) {
        this->errtypes.push_back(std::vector<uint8_t>());
        this->sync_group.push_back(std::vector<int>());
        this->callq.push_back(std::vector<float>());
        this->credit.push_back(std::vector<float>());
        this->ref_ed.push_back(std::vector<int>());
        this->query_ed.push_back(std::vector<int>());
    }
}


/******************************************************************************/


/* Efficiently remove many elements from a vector. */
template<typename T>
void _remove_indices(std::vector<T> & vector, const std::vector<int> & to_remove)
{
    auto vector_base = vector.begin();
    size_t down_by = 0;

    for (auto iter = to_remove.cbegin(); iter < to_remove.cend(); iter++, down_by++)
        {
        size_t next = (iter + 1 == to_remove.cend() ? vector.size() : *(iter + 1));
        std::move(vector_base + *iter + 1, vector_base + next, vector_base + *iter - down_by);
    }
    vector.resize(vector.size() - to_remove.size());
}


/* Remove many variants at once. */
void ctgVariants::remove_vars(const std::vector<int> & indices) {
    // set for all variants
    _remove_indices(this->poss, indices);
    _remove_indices(this->rlens, indices);
    _remove_indices(this->types, indices);
    _remove_indices(this->locs, indices);
    _remove_indices(this->refs, indices);
    _remove_indices(this->alts, indices);
    _remove_indices(this->orig_gts, indices);
    _remove_indices(this->gt_quals, indices);
    _remove_indices(this->var_quals, indices);
    _remove_indices(this->phase_sets, indices);
    this->n -= int(indices.size());

    // added during precision/recall analysis
    _remove_indices(this->calc_gts, indices);
    for (int hi = 0; hi < HAPS; hi++) {
        _remove_indices(this->errtypes[hi], indices);
        _remove_indices(this->sync_group[hi], indices);
        _remove_indices(this->callq[hi], indices);
        _remove_indices(this->ref_ed[hi], indices);
        _remove_indices(this->query_ed[hi], indices);
        _remove_indices(this->credit[hi], indices);
    }

    // added during phasing analysis
    _remove_indices(this->phases, indices);
    _remove_indices(this->pb_phases, indices);
    _remove_indices(this->ac_errtype, indices);
}

void ctgVariants::add_var(std::shared_ptr<ctgVariants> other_vars, int idx) {
    this->add_var(
        other_vars->poss[idx], 
        other_vars->rlens[idx], 
        other_vars->types[idx], 
        other_vars->locs[idx],
        other_vars->refs[idx],
        other_vars->alts[idx],
        other_vars->orig_gts[idx],
        other_vars->gt_quals[idx],
        other_vars->var_quals[idx],
        other_vars->phase_sets[idx]);
}

void ctgVariants::add_var(int pos, int rlen, uint8_t type, uint8_t loc,
        std::string ref, std::string alt, uint8_t orig_gt, float gq, float vq, int ps) {
    // set for all variants
    this->poss.push_back(pos);
    this->rlens.push_back(rlen);
    this->types.push_back(type);
    this->locs.push_back(loc);
    this->refs.push_back(ref);
    this->alts.push_back(alt);
    this->orig_gts.push_back(orig_gt);
    this->gt_quals.push_back(gq);
    this->var_quals.push_back(std::min(vq, float(g.max_qual)));
    this->phase_sets.push_back(ps);
    this->n++;

    // added during precision/recall analysis
    this->calc_gts.push_back(GT_REF_REF);
    for (int hi = 0; hi < HAPS; hi++) {
        this->errtypes[hi].push_back(ERRTYPE_UN);
        this->sync_group[hi].push_back(0);
        this->callq[hi].push_back(0);
        this->ref_ed[hi].push_back(0);
        this->query_ed[hi].push_back(0);
        this->credit[hi].push_back(0);
    }

    // added during phasing analysis
    this->phases.push_back(PHASE_NONE);
    this->pb_phases.push_back(PHASE_NONE);
    this->ac_errtype.push_back(AC_UNKNOWN);
}


/******************************************************************************/


int ctgVariants::get_vartype(int vi) { 
    if (this->types[vi] == TYPE_SUB) { // SNP
        return VARTYPE_SNP;
    } else if ((this->types[vi] == TYPE_INS && // small INDEL
                int(this->alts[vi].size()) < g.sv_threshold) ||
            (this->types[vi] == TYPE_DEL &&
             int(this->refs[vi].size()) < g.sv_threshold)) {
        return VARTYPE_INDEL;
    } else { // SV
        return VARTYPE_SV;
    }
}


/******************************************************************************/


/* Precision/Recall calculation allows dynamically adjusting the variant genotype.
 * Afterwards, this functions records how the genotype allele count has changed.
 */
int ctgVariants::set_allele_errtype(int vi) { 
    if (this->calc_gts[vi] == GT_ALT1_REF || this->calc_gts[vi] == GT_REF_ALT1) {
        if (this->orig_gts[vi] == GT_ALT1_ALT1) {
            return this->ac_errtype[vi] = AC_ERR_1_TO_2;
        } else if (this->orig_gts[vi] == GT_ALT1_REF || this->orig_gts[vi] == GT_REF_ALT1) {
            return this->ac_errtype[vi] = AC_ERR_1_TO_1;
        } else {
            return this->ac_errtype[vi] = AC_ERR_1_TO_0;
        }
    } else if (this->calc_gts[vi] == GT_ALT1_ALT1) {
        if (this->orig_gts[vi] == GT_ALT1_ALT1) {
            return this->ac_errtype[vi] = AC_ERR_2_TO_2;
        } else if (this->orig_gts[vi] == GT_ALT1_REF || this->orig_gts[vi] == GT_REF_ALT1) {
            return this->ac_errtype[vi] = AC_ERR_2_TO_1;
        } else {
            return this->ac_errtype[vi] = AC_ERR_2_TO_0;
        }
    } else if (this->calc_gts[vi] == GT_REF_REF) {
        if (this->orig_gts[vi] == GT_ALT1_ALT1) {
            return this->ac_errtype[vi] = AC_ERR_0_TO_2;
        } else if (this->orig_gts[vi] == GT_ALT1_REF || this->orig_gts[vi] == GT_REF_ALT1) {
            return this->ac_errtype[vi] = AC_ERR_0_TO_1;
        }
    }
    return AC_UNKNOWN;
}


/******************************************************************************/


void variantData::print_phase_info(int callset) {

    int total_phase_sets = 0;
    std::vector<int> phase_set_sizes;
    int ctg_idx = 0;

    for (const std::string & ctg : this->contigs) { // for each contig

        // get first phase set (to backfill all preceding zeros)
        int first_pos = std::numeric_limits<int>::max();
        int first_phase_set = 0;
        for (int hi = 0; hi < HAPS; hi++) {
            std::shared_ptr<ctgVariants> vars = this->variants[hi][ctg];
            for (int vi = 0; vi < vars->n; vi++) {
                if (vars->phase_sets[vi] != 0) {
                    if (vars->poss[vi] < first_pos) {
                        first_pos = vars->poss[vi]; 
                        first_phase_set = vars->phase_sets[vi];
                    }
                    break;
                }
            }
        }

        // exit early if this contig has no phase sets
        if (first_phase_set == 0) {
            phase_set_sizes.push_back(this->lengths[ctg_idx]);
            total_phase_sets++;
            continue;
        } else { // set phase set up until first PS
            for (int hi = 0; hi < HAPS; hi++) {
                std::shared_ptr<ctgVariants> vars = this->variants[hi][ctg];
                for (int vi = 0; vi < vars->n; vi++) {
                    if (vars->phase_sets[vi] != 0) {
                        break;
                    }
                    vars->phase_sets[vi] = first_phase_set;
                }
            }
        }

        // TODO: why is hi needed here? when are variants merged? shouldn't this iterate over one list?
        std::vector<int> vi(HAPS, 0);
        int ps_beg = 0; int ps_end = 0;
        int phase_set = 0;
        while (vi[0] < this->variants[0][ctg]->n || vi[1] < this->variants[1][ctg]->n) {

            // select haplotype with variant that occurs first
            int hi = 0;
            if (vi[0] >= this->variants[0][ctg]->n) {
                hi = 1;
            } else if (vi[1] >= this->variants[1][ctg]->n) {
                hi = 0;
            } else {
                hi = this->variants[0][ctg]->poss[vi[0]] > this->variants[1][ctg]->poss[vi[1]];
            }
            
            // TODO: for now, only start new phase set if new PS tag found
            std::shared_ptr<ctgVariants> vars = this->variants[hi][ctg];
            if (vars->phase_sets[vi[hi]] != 0) { // variant is phased
                if (vars->phase_sets[vi[hi]] != phase_set) { // new phase set, save old
                    if (phase_set) phase_set_sizes.push_back(ps_end - ps_beg);
                    phase_set = vars->phase_sets[vi[hi]];
                    ps_beg = vars->poss[vi[hi]];
                    ps_end = vars->poss[vi[hi]] + vars->rlens[vi[hi]];
                    total_phase_sets++;
                } else { // same 
                    ps_end = std::max(ps_end, vars->poss[vi[hi]] + vars->rlens[vi[hi]]);
                }
            } else { // set unphased variant phase set to current phase set
                // TODO: only set phase sets for 1|1 variants (for when we add unphased eval)
                vars->phase_sets[vi[hi]] = phase_set;
            }
            
            vi[hi]++;
        }

        // add final phase set on contig
        phase_set_sizes.push_back(ps_end - ps_beg);
        total_phase_sets++;
        ctg_idx++;
    }

    // calculate phaseset NG50
    size_t total_bases = 0;
    for (size_t i = 0; i < this->contigs.size(); i++) {
        total_bases += lengths[i];
    }

    int pb_ng50 = calc_ng50(phase_set_sizes, total_bases);

    if (g.verbosity >= 1) INFO("          %s phase sets: %d",
            callset_strs[callset].data(), total_phase_sets);
    if (g.verbosity >= 1) INFO("    %s phase block NG50: %d", 
            callset_strs[callset].data(), pb_ng50);
    if (g.verbosity >= 1) INFO("        Total contig bases: %zu", total_bases);
}


/*******************************************************************************/


void variantData::write_vcf(std::string out_vcf_fn) {

    // VCF header
    FILE* out_vcf = fopen(out_vcf_fn.data(), "w");
    const std::chrono::time_point now{std::chrono::system_clock::now()};
    time_t tt = std::chrono::system_clock::to_time_t(now);
    tm local_time = *localtime(&tt);
    fprintf(out_vcf, "##fileformat=VCFv4.2\n");
    fprintf(out_vcf, "##fileDate=%04d%02d%02d\n", local_time.tm_year + 1900, 
            local_time.tm_mon + 1, local_time.tm_mday);
    for (size_t i = 0; i < this->contigs.size(); i++)
        fprintf(out_vcf, "##contig=<ID=%s,length=%d,ploidy=%d>\n", 
                this->contigs[i].data(), this->lengths[i], this->ploidy[i]);
    fprintf(out_vcf, "##FILTER=<ID=PASS,Description=\"All filters passed\">\n");
    fprintf(out_vcf, "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n");
    fprintf(out_vcf, "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t%s\n",
            this->sample.data());

    // write variants
    for (std::string ctg : this->contigs) {
        std::vector<size_t> ptrs = {0, 0};
        int p = this->ploidy[std::find(contigs.begin(), contigs.end(), ctg) - contigs.begin()];
        while (ptrs[HAP1] < this->variants[HAP1][ctg]->poss.size() ||
                ptrs[HAP2] < this->variants[HAP2][ctg]->poss.size()) {

            // get next positions, set flags for which haps
            int pos_hap1 = ptrs[HAP1] < this->variants[HAP1][ctg]->poss.size() ? 
                this->variants[HAP1][ctg]->poss[ptrs[HAP1]] : std::numeric_limits<int>::max();
            int pos_hap2 = ptrs[HAP2] < this->variants[HAP2][ctg]->poss.size() ? 
                this->variants[HAP2][ctg]->poss[ptrs[HAP2]] : std::numeric_limits<int>::max();

            // indels include previous base, adjust position
            if (ptrs[HAP1] < this->variants[HAP1][ctg]->types.size() && 
                    (this->variants[HAP1][ctg]->types[ptrs[HAP1]] == TYPE_INS || 
                    this->variants[HAP1][ctg]->types[ptrs[HAP1]] == TYPE_DEL)) pos_hap1--;
            if (ptrs[HAP2] < this->variants[HAP2][ctg]->types.size() && 
                    (this->variants[HAP2][ctg]->types[ptrs[HAP2]] == TYPE_INS || 
                    this->variants[HAP2][ctg]->types[ptrs[HAP2]] == TYPE_DEL)) pos_hap2--;
            int pos = std::min(pos_hap1, pos_hap2);
            bool hap1 = (pos_hap1 == pos);
            bool hap2 = (pos_hap2 == pos);

            // add variants to output VCF file
            if (hap1 && hap2) {
                if (this->variants[HAP1][ctg]->refs[ptrs[HAP1]] == 
                        this->variants[HAP2][ctg]->refs[ptrs[HAP2]] &&
                        this->variants[HAP1][ctg]->alts[ptrs[HAP1]] == 
                        this->variants[HAP2][ctg]->alts[ptrs[HAP2]]) {
                    
                    // homozygous variant (1|1)
                    print_variant(out_vcf, ctg, pos, 
                            this->variants[HAP1][ctg]->types[ptrs[HAP1]],
                            this->variants[HAP1][ctg]->refs[ptrs[HAP1]],
                            this->variants[HAP1][ctg]->alts[ptrs[HAP1]],
                            this->variants[HAP1][ctg]->var_quals[ptrs[HAP1]], "1|1");
                    
                } else {
                    // two separate phased variants (0|1 + 1|0)
                    print_variant(out_vcf, ctg, pos, 
                            this->variants[HAP1][ctg]->types[ptrs[HAP1]],
                            this->variants[HAP1][ctg]->refs[ptrs[HAP1]],
                            this->variants[HAP1][ctg]->alts[ptrs[HAP1]],
                            this->variants[HAP1][ctg]->var_quals[ptrs[HAP1]], "1|0");
                    print_variant(out_vcf, ctg, pos, 
                            this->variants[HAP2][ctg]->types[ptrs[HAP2]],
                            this->variants[HAP2][ctg]->refs[ptrs[HAP2]],
                            this->variants[HAP2][ctg]->alts[ptrs[HAP2]],
                            this->variants[HAP2][ctg]->var_quals[ptrs[HAP2]], "0|1");
                }

            } else if (hap1) { // 1|0
                print_variant(out_vcf, ctg, pos, 
                        this->variants[HAP1][ctg]->types[ptrs[HAP1]],
                        this->variants[HAP1][ctg]->refs[ptrs[HAP1]],
                        this->variants[HAP1][ctg]->alts[ptrs[HAP1]],
                        this->variants[HAP1][ctg]->var_quals[ptrs[HAP1]], p == 1 ? "1" : "1|0");

            } else if (hap2) { // 0|1
                print_variant(out_vcf, ctg, pos, 
                        this->variants[HAP2][ctg]->types[ptrs[HAP2]],
                        this->variants[HAP2][ctg]->refs[ptrs[HAP2]],
                        this->variants[HAP2][ctg]->alts[ptrs[HAP2]],
                        this->variants[HAP2][ctg]->var_quals[ptrs[HAP2]],  p == 1 ? "1" :"0|1");
            }

            // update pointers
            if (hap1) ptrs[HAP1]++;
            if (hap2) ptrs[HAP2]++;
        }
    }
    fclose(out_vcf);
}


/*******************************************************************************/
/* This function is used to map haplotypes between orig_gt and calc_gt, in order to determine which
   calc_gt data should be reported for the initial orig_gt variant, and is complicated by the
   fact that the two may contain a different number of non-reference alleles. 
 */ 
bool ctgVariants::calcgt_is_swapped(int vi /* variant index */) const {
    // 0|0,0|0 and 0|1,0|1 and 1|0,1|0 and 1|1,1|1
    if (this->orig_gts[vi] == this->calc_gts[vi]) {
        return false;
    }
    // orig_gt == 1|1 or 0|0 or calc_gt == 0|0, all/no calc_gt data will be reported, order doesn't matter
    else if (this->orig_gts[vi] == GT_ALT1_ALT1 || 
             this->orig_gts[vi] == GT_REF_REF || 
             this->calc_gts[vi] == GT_REF_REF) {
        return false;
    }
    // 0|1,1|0 and 1|0,0|1
    else if ((this->orig_gts[vi] == GT_REF_ALT1 && this->calc_gts[vi] == GT_ALT1_REF) || 
             (this->orig_gts[vi] == GT_ALT1_REF && this->calc_gts[vi] == GT_REF_ALT1)) {
        return true;
    }
    // orig_gt = 0|1, choose better calc_gt
    else if (this->orig_gts[vi] == GT_REF_ALT1 && this->calc_gts[vi] == GT_ALT1_ALT1) {
        return this->credit[HAP1][vi] > this->credit[HAP2][vi];
    }
    // orig_gt = 1|0, choose better calc_gt
    else if (this->orig_gts[vi] == GT_ALT1_REF && this->calc_gts[vi] == GT_ALT1_ALT1) {
        return this->credit[HAP2][vi] > this->credit[HAP1][vi];
    } else {
        ERROR("Unexpected orig/calc genotypes for variant (%s -> %s) at pos %d: orig=%s calc=%s",
                this->refs[vi].data(),
                this->alts[vi].data(),
                this->poss[vi],
                gt_strs[this->orig_gts[vi]].data(),
                gt_strs[this->calc_gts[vi]].data()
        );
    }
}

bool ctgVariants::var_on_hap(int var_idx, int hap, bool calc) const {
    int gt = calc ? this->calc_gts[var_idx] : this->orig_gts[var_idx]; // simple gt, always (0|1, 1|0, or 1|1)
    if (hap == 0 && (gt == GT_ALT1 || gt == GT_ALT1_REF || gt == GT_ALT1_ALT1))
        return true;
    if (hap == 1 && (gt == GT_ALT1 || gt == GT_REF_ALT1 || gt == GT_ALT1_ALT1))
        return true;
    if (hap > 1)
        ERROR("Unexpected haplotype %d for variant (%s -> %s) at pos %d", int(hap), 
                this->refs[var_idx].data(),
                this->alts[var_idx].data(),
                this->poss[var_idx]);
    return false;
}

/*******************************************************************************/

/* Set (or unset) the presence of an alternate allele in a variant's genotype.
 */
void ctgVariants::set_var_calcgt_on_hap(int var_idx, int hap, bool set, bool ignore_errors) {
    if (hap > 1) ERROR("Unexpected hap idx %d in set_var_calcgt_on_hap()", hap);

    if (this->calc_gts[var_idx] == GT_REF_REF) {
        if (set) {
            this->calc_gts[var_idx] = hap == 0 ? GT_ALT1_REF : GT_REF_ALT1;
        } else { // unset
            if (!ignore_errors) ERROR("Variant calc_gt already unset for variant %d at %s:%d hap %d",
                    var_idx, this->ctg.data(), this->poss[var_idx], hap);
        }

    } else if (this->calc_gts[var_idx] == GT_REF_ALT1) {
        if (set) {
            if (hap == 1) {
                if (!ignore_errors) ERROR("Variant calc_gt already set for variant %d at %s:%d hap %d",
                    var_idx, this->ctg.data(), this->poss[var_idx], hap);
            } else {
                this->calc_gts[var_idx] = GT_ALT1_ALT1;
            }
        } else { // unset
            if (hap == 0) {
                if (!ignore_errors) ERROR("Variant calc_gt already unset for variant %d at %s:%d hap %d",
                    var_idx, this->ctg.data(), this->poss[var_idx], hap);
            } else {
                this->calc_gts[var_idx] = GT_REF_REF;
            }
        }

    } else if (this->calc_gts[var_idx] == GT_ALT1_REF) {
        if (set) {
            if (hap == 0) {
                if (!ignore_errors) ERROR("Variant calc_gt already set for variant %d at %s:%d hap %d",
                    var_idx, this->ctg.data(), this->poss[var_idx], hap);
            } else {
                this->calc_gts[var_idx] = GT_ALT1_ALT1;
            }
        } else { // unset
            if (hap == 1) {
                if (!ignore_errors) ERROR("Variant calc_gt already unset for variant %d at %s:%d hap %d",
                    var_idx, this->ctg.data(), this->poss[var_idx], hap);
            } else {
                this->calc_gts[var_idx] = GT_REF_REF;
            }
        }

    } else if (this->calc_gts[var_idx] == GT_ALT1_ALT1) {
        if (set) {
            if (!ignore_errors) ERROR("Variant calc_gt already set for variant %d at %s:%d hap %d",
                    var_idx, this->ctg.data(), this->poss[var_idx], hap);
        } else {
            this->calc_gts[var_idx] = hap == 0 ? GT_REF_ALT1 : GT_ALT1_REF;
        }

    } else {
        ERROR("Unexpected calc_gts value '%s' in set_var_calcgt_on_hap() for variant %d at %s:%d", 
            gt_strs[this->calc_gts[var_idx]].data(), var_idx, this->ctg.data(), this->poss[var_idx]);
    }
}

/*******************************************************************************/

void ctgVariants::print_var_info(FILE* out_fp, std::shared_ptr<fastaData> ref, 
        const std::string & ctg, int idx) {
    char ref_base;
    switch (this->types[idx]) {
    case TYPE_SUB:
        fprintf(out_fp, "%s\t%d\t.\t%s\t%s\t.\tPASS\t.\tGT:BD:BC:RD:QD:BK:QQ:SC:SG:PS:PB:BS:VP:FE:GE", 
                ctg.data(), this->poss[idx]+1, this->refs[idx].data(), 
                this->alts[idx].data());
        break;
    case TYPE_INS:
    case TYPE_DEL:
        ref_base = ref->fasta.at(ctg)[this->poss[idx]-1];
        fprintf(out_fp, "%s\t%d\t.\t%s\t%s\t.\tPASS\t.\tGT:BD:BC:RD:QD:BK:QQ:SC:SG:PS:PB:BS:VP:FE:GE", ctg.data(), 
                this->poss[idx], (ref_base + this->refs[idx]).data(), 
                (ref_base + this->alts[idx]).data());
        break;
    default:
        ERROR("print_variant not implemented for type %d", this->types[idx]);
    }
}


void ctgVariants::print_var_empty(FILE* out_fp, int sc_idx, 
        int phase_block, bool query /* = false */) {
    fprintf(out_fp, "\t.:.:.:.:.:.:.:%d:.:.:%d:.:.:.:.%s", sc_idx, phase_block, query ? "\n" : "");
}


void ctgVariants::print_var_sample(FILE* out_fp, int vi, int hi, const std::string & gt,
        int sc_idx, int phase_block, bool phase_switch, bool phase_flip, bool query /* = false */) {

    // get categorization
    std::string errtype;
    std::string match_type;
    if (this->credit[hi][vi] == 1) {
        errtype = "TP"; match_type = "gm";
    } else if (this->credit[hi][vi] == 0) {
        errtype = query ? "FP" : "FN"; match_type = ".";
    } else if (this->credit[hi][vi] >= g.credit_threshold) {
        errtype = "TP"; match_type = "lm";
    } else {
        errtype = query ? "FP" : "FN"; match_type = "lm";
    }

    fprintf(out_fp, "\t%s:%s:%f:%s:%s:%s:%d:%d:%d:%d:%d:%s:%s:%s:%s%s", gt.data(), errtype.data(), 
            this->credit[hi][vi], 
            this->ref_ed[hi][vi] == 0 ? "." : 
                std::to_string(this->ref_ed[hi][vi]).data(),
            this->ref_ed[hi][vi] == 0 ? "." : 
                std::to_string(this->query_ed[hi][vi]).data(),
            match_type.data(), int(this->var_quals[vi]), sc_idx, 
            int(this->sync_group[hi][vi]), this->phase_sets[vi], phase_block,
            query ? (phase_switch ? "1" : "0") : "." , 
            phase_strs[this->phases[vi]].data(),
            query ? (phase_flip ? "1" : "0") : "." , 
            ac_strs[this->ac_errtype[vi]].data(),
            query ? "\n" : "");
}


/*******************************************************************************/


void variantData::print_variant(FILE* out_fp, const std::string & ctg, int pos, int type,
        const std::string & ref, const std::string & alt, float qual, const std::string & gt) {

    char ref_base;
    switch (type) {
    case TYPE_SUB:
        fprintf(out_fp, "%s\t%d\t.\t%s\t%s\t%f\tPASS\t.\tGT\t%s\n", ctg.data(),
            pos+1, ref.data(), alt.data(), qual, gt.data());
        break;
    case TYPE_INS:
    case TYPE_DEL:
        try {
            ref_base = this->ref->fasta.at(ctg)[pos];
            fprintf(out_fp, "%s\t%d\t.\t%s\t%s\t%f\tPASS\t.\tGT\t%s\n", ctg.data(), 
                    pos+1, (ref_base + ref).data(), (ref_base + alt).data(), 
                    qual, gt.data());
        } catch (const std::out_of_range & e) {
            ERROR("Contig '%s' not in reference FASTA (print_variant)", ctg.data());
        }
        break;
    default:
        ERROR("print_variant not implemented for type %d", type);
    }
}

void variantData::set_header(const std::shared_ptr<variantData> vcf) {
    this->filename = vcf->filename;
    this->sample = vcf->sample;
    this->contigs = vcf->contigs;
    this->lengths = vcf->lengths;
    this->ploidy = vcf->ploidy;
    this->ref = vcf->ref;
    for (const std::string & ctg : this->contigs)
        for (int hap = 0; hap < 2; hap++)
            this->variants[hap][ctg] = 
                std::shared_ptr<ctgVariants>(new ctgVariants(ctg));
}


/* Copy variants from `cigar` string to `variantData`. */
void variantData::add_variants(
        const std::vector<int> & cigar, 
        int hap, int ref_pos,
        const std::string & ctg, 
        const std::string & query, 
        const std::string & ref, 
        int qual, int phase_set) {

    int query_idx = 0;
    int ref_idx = 0;
    for (size_t cig_idx = 0; cig_idx < cigar.size(); ) {
        int indel_len = 0;
        switch (cigar[cig_idx]) {

            case PTR_MAT: // no variant, update pointers
                cig_idx += 2;
                ref_idx++;
                query_idx++;
                break;

            case PTR_SUB: // substitution
                cig_idx += 2;
                this->variants[hap][ctg]->add_var(ref_pos+ref_idx, 1, 
                        TYPE_SUB, BED_INSIDE, std::string(1,ref[ref_idx]), 
                        std::string(1,query[query_idx]), 
                        GT_REF_REF, g.max_qual, qual, phase_set);
                ref_idx++;
                query_idx++;
                break;

            case PTR_DEL: // deletion
                cig_idx++; indel_len++;

                // multi-base deletion
                while (cig_idx < cigar.size() && cigar[cig_idx] == PTR_DEL) {
                    cig_idx++; indel_len++;
                }
                this->variants[hap][ctg]->add_var(ref_pos+ref_idx,
                        indel_len, TYPE_DEL, BED_INSIDE,
                        ref.substr(ref_idx, indel_len),
                        "", GT_REF_REF, g.max_qual, qual, phase_set);
                ref_idx += indel_len;
                break;

            case PTR_INS: // insertion
                cig_idx++; indel_len++;

                // multi-base insertion
                while (cig_idx < cigar.size() && cigar[cig_idx] == PTR_INS) {
                    cig_idx++; indel_len++;
                }
                this->variants[hap][ctg]->add_var(ref_pos+ref_idx,
                        0, TYPE_INS, BED_INSIDE, "", 
                        query.substr(query_idx, indel_len), 
                        GT_REF_REF, g.max_qual, qual, phase_set);
                query_idx += indel_len;
                break;
        }
    }
}

/******************************************************************************/

variantData::variantData() : callset(QUERY), variants(HAPS) { ; }

void parse_variants(const std::string & vcf_fn, 
        std::shared_ptr<variantData> variant_data,
        std::shared_ptr<variantData> large_variant_data,
        std::shared_ptr<fastaData> reference, 
        int callset) {

    // set reference fasta pointer
    variant_data->ref = reference;
    large_variant_data->ref = reference;
    variant_data->filename = vcf_fn;
    large_variant_data->filename = vcf_fn;

    if (callset < 0 || callset >= CALLSETS)
        ERROR("Invalid callset (%d).", callset);
    variant_data->callset = callset;
    large_variant_data->callset = callset;

    if (g.verbosity >= 1) INFO(" ");
    if (g.verbosity >= 1) INFO("%s[%s %d/%d] Parsing %s VCF%s '%s'", COLOR_PURPLE,
            callset == QUERY ? "Q" : "T", TIME_READ, TIME_TOTAL-1, callset_strs[callset].data(), 
            COLOR_WHITE, vcf_fn.data());
    htsFile* vcf = bcf_open(vcf_fn.data(), "r");

    // counters
    int nctg   = 0;                     // number of ctgs
    std::vector< std::vector<int> > 
        ntypes(2, std::vector<int>(type_strs.size(), 0));
    int n      = 0;                     // total number of records in file
    std::vector<int> npass  = {0, 0};   // records PASSing all filters

    // data
    bool print = g.verbosity >= 1;
    std::vector<int> prev_end = {-g.cluster_min_gap*2, -g.cluster_min_gap*2};
    std::vector<int> prev_type = {TYPE_SUB, TYPE_SUB};
    std::unordered_map<int, bool> prev_rids;
    int prev_rid = -1;
    std::unordered_map<int, int> ctglens;
    std::string ctg;
    std::vector<int> nregions(region_strs.size(), 0);
    std::vector<int> pass_min_qual = {FALSE, FALSE};

    // quality data for each call
    int GQ_memsize = 0;
    int ngq     = 0;
    int * gq    = (int*) malloc(sizeof(int));
    float * fgq = (float*) malloc(sizeof(float));
    bool int_qual = true;

    // genotype data for each call
    int GT_memsize   = 0;
    int ngt       = 0;
    std::vector<int> GT_counts(gt_strs.size(), 0);
    int * gt      = NULL;
    bool gt_warn  = false;

    // phase set data for each call
    int PS_memsize   = 0;
    int nPS       = 0;
    int * PS      = NULL;
    bool PS_warn  = false;
    int phase_set = 0;

    /* int gq_missing_total = 0; */
    int PS_missing_total = 0;
    int overlapping_var_total = 0;
    int spanning_del_total = 0;
    int unknown_allele_total = 0;
    int unphased_gt_total = 0;
    int large_var_total = 0;
    int too_large_var_total = 0;
    int wrong_ploidy_total = 0;
    int multi_total = 0;
    int ref_call_total = 0;
    int complex_total = 0;
    int failed_filter_total = 0;
    
    // read header
    bcf1_t * rec  = NULL;
    bcf_hdr_t *hdr = bcf_hdr_read(vcf);
    bool pass = false;
    for (int i = 0; i < hdr->nhrec; i++) { // for each header record (line)

        /* // DEBUG HEADER PRINT */
        /* printf("%s=%s\n", hdr->hrec[i]->key, hdr->hrec[i]->value); */
        /* for (int j = 0; j < hdr->hrec[i]->nkeys; j++) { */
        /*     printf("  %s=%s\n", hdr->hrec[i]->keys[j], hdr->hrec[i]->vals[j]); */
        /* } */

        // search all FILTER lines
        if (hdr->hrec[i]->type == BCF_HL_FLT) {

            // is this a filter we selected?
            int filter_idx = -1;
            for (int j = 0; j < hdr->hrec[i]->nkeys; j++) { // for each key in this FILTER record
                for (int fi = 0; fi < int(g.filters.size()); fi++) { // for each filter
                    if (std::string(hdr->hrec[i]->keys[j]) == std::string("ID") &&  
                            std::string(hdr->hrec[i]->vals[j]) == g.filters[fi]) {
                        filter_idx = fi;
                    }
                }
            }

            // if so, save this filter's index (to compare against for variants)
            if (filter_idx >= 0) {
                for (int j = 0; j < hdr->hrec[i]->nkeys; j++) {
                    if (std::string(hdr->hrec[i]->keys[j]) == std::string("IDX")) {
                        g.filter_ids[filter_idx] = std::stoi(hdr->hrec[i]->vals[j]);
                    }
                }
            }
        }

        // store contig lengths (for output VCF)
        else if (hdr->hrec[i]->type == BCF_HL_CTG) {
            int length = -1;
            int idx = -1;
            for (int j = 0; j < hdr->hrec[i]->nkeys; j++) {
                if (hdr->hrec[i]->keys[j] == std::string("IDX"))
                    idx = std::stoi(hdr->hrec[i]->vals[j]);
                else if (hdr->hrec[i]->keys[j] == std::string("length"))
                    length = std::stoi(hdr->hrec[i]->vals[j]);
            }
            if (length >= 0 && idx >= 0) {
                ctglens[idx] = length;
            } else {
                ERROR("%s VCF header contig line didn't have 'IDX' and 'length'",
                        callset_strs[callset].data());
            }
        }
    }

    if (bcf_hdr_nsamples(hdr) != 1) 
        ERROR("Expected 1 sample but found %d in %s VCF '%s'", bcf_hdr_nsamples(hdr),
                callset_strs[callset].data(), vcf_fn.data());
    variant_data->sample = hdr->samples[0];
    large_variant_data->sample = hdr->samples[0];

    // verify that the filters we selected were in the VCF
    for (int fi = 0; fi < int(g.filters.size()); fi++) {
        if (g.filter_ids[fi] < 0)
            WARN("Filter '%s' not found in %s VCF", g.filters[fi].data(),
                    callset_strs[callset].data());
    }

    // report names of all the ctgs in the VCF file
    const char **ctgnames = NULL;
    ctgnames = bcf_hdr_seqnames(hdr, &nctg);
    if (ctgnames == NULL) {
        ERROR("Failed to read %s VCF '%s' header", 
                callset_strs[callset].data(), vcf_fn.data());
        goto error1;
    }
    for(int i = 0; i < nctg; i++) {
        variant_data->variants[HAP1][ctgnames[i]] = 
                std::shared_ptr<ctgVariants>(new ctgVariants(ctgnames[i]));
        variant_data->variants[HAP2][ctgnames[i]] = 
                std::shared_ptr<ctgVariants>(new ctgVariants(ctgnames[i]));
        large_variant_data->variants[HAP1][ctgnames[i]] = 
                std::shared_ptr<ctgVariants>(new ctgVariants(ctgnames[i]));
        large_variant_data->variants[HAP2][ctgnames[i]] = 
                std::shared_ptr<ctgVariants>(new ctgVariants(ctgnames[i]));
    }

    // struct for storing each record
    rec = bcf_init();
    if (rec == NULL) {
        ERROR("Failed to read %s VCF '%s' records", 
                callset_strs[callset].data(), vcf_fn.data());
        goto error2;
    }
    
    while (bcf_read(vcf, hdr, rec) == 0) {

        ctg = ctgnames[rec->rid];
        if (rec->rid != prev_rid) {

            // start new contig
            prev_rid = rec->rid;
            if (prev_rids.find(rec->rid) != prev_rids.end()) {
                ERROR("Unsorted %s VCF '%s', contig '%s' already parsed", 
                        callset_strs[callset].data(), vcf_fn.data(), ctg.data());
            } else {
                variant_data->contigs.push_back(ctg);
                variant_data->ploidy.push_back(0);
                variant_data->lengths.push_back(ctglens[rec->rid]);
                large_variant_data->contigs.push_back(ctg);
                large_variant_data->ploidy.push_back(0);
                large_variant_data->lengths.push_back(ctglens[rec->rid]);
                prev_end = {-g.cluster_min_gap*2, -g.cluster_min_gap*2};
                prev_type = {TYPE_SUB, TYPE_SUB};
            }
        }

        // unpack info (populates rec->d allele info)
        bcf_unpack(rec, BCF_UN_ALL);
        n++;

        // check that variant contained a passing filter
        pass = false;
        if (rec->d.n_flt == 0 || g.filters.size() == 0) { // no filters, default pass
            pass = true;
        } else {
            for (int i = 0; i < rec->d.n_flt; i++) { // check this variant's filters
                for (int fi = 0; fi < int(g.filters.size()); fi++) {
                    if (rec->d.flt[i] == g.filter_ids[fi]) pass = true;
                }
            }
        }

        // variant doesn't contain a passing filter
        if (!pass) {
            failed_filter_total++;
            continue;
        }

        // check that variant exceeds min_qual
        float vq = rec->qual;
        if (std::isnan(vq)) vq = 0; // no quality reported (.)
        pass = vq >= g.min_qual;
        pass_min_qual[pass]++;
        if (!pass) continue;

        // parse GQ in either INT or FLOAT format
        if (int_qual) {
            ngq = bcf_get_format_int32(hdr, rec, "GQ", &gq, &GQ_memsize);
            if (ngq == -2) {
                ngq = bcf_get_format_float(hdr, rec, "GQ", &fgq, &GQ_memsize);
                gq[0] = int(fgq[0]);
                int_qual = false;
            }
        }
        else {
            ngq = bcf_get_format_float(hdr, rec, "GQ", &fgq, &GQ_memsize);
            gq[0] = int(fgq[0]);
        }
        if ( ngq == -3 || ngq == -1 ) { // missing
            /* if (g.verbosity > 1 || !gq_missing_total) */
            /*     WARN("No GQ tag in %s VCF at %s:%lld", */
            /*             callset_strs[callset].data(), ctg.data(), (long long)rec->pos); */
            /* gq_missing_total++; // only warn once */
            gq[0] = 0;
        } else if (ngq < 0) { // other error
            ERROR("Failed to read %s GQ at %s:%lld", 
                    callset_strs[callset].data(), ctg.data(), (long long)rec->pos);
        }

        // parse GT: https://github.com/samtools/htslib/blob/99415e2a2ce26bdbf4e910954330ea769de2c3f0/htslib/vcf.h#L1096
        ngt = bcf_get_format_int32(hdr, rec, "GT", &gt, &GT_memsize);
        if (ngt == -1) { // GT not defined in header
            if (!gt_warn) {
                gt_warn = true;
                WARN("'GT' tag not defined in %s VCF header, assuming monoploid", 
                        callset_strs[callset].data());
            }
        } else if (ngt <= 0) { // other error
            ERROR("Failed to read %s GT at %s:%lld", 
                    callset_strs[callset].data(), ctg.data(), (long long)rec->pos);
        }

        // update ploidy info
        int ctg_idx = std::find(variant_data->contigs.begin(), variant_data->contigs.end(), ctg) 
                - variant_data->contigs.begin();
        if (variant_data->ploidy[ctg_idx] != 0) { // already set, enforce it doesn't change
            if (std::abs(ngt) != variant_data->ploidy[ctg_idx] && ctg[ctg.size()-1] != 'X') {
                if (g.verbosity > 1)
                    WARN("Expected ploidy %d for all variants on contig '%s',"
                          " found ploidy %d at %s:%lld in %s VCF.", variant_data->ploidy[ctg_idx],
                        ctg.data(), std::abs(ngt), ctg.data(), (long long)rec->pos,
                        callset_strs[callset].data());
                wrong_ploidy_total += 1;
            }
        } else { // set ploidy for this contig
            variant_data->ploidy[ctg_idx] = std::abs(ngt);
            large_variant_data->ploidy[ctg_idx] = std::abs(ngt);
        }

        // parse genotype info
        int orig_gt = GT_REF_REF;
        bool same = false;
        if (ngt == -1) { // no info, assume monoploid
            orig_gt = GT_ALT1;
        } else if (ngt == 1) { // monoploid/haploid

            // set 1 if allele_idx > 0
            orig_gt = bcf_gt_allele(gt[0]) ? GT_ALT1 : GT_REF;

        } else if (ngt == 2) { // diploid

            // missing, ignore
            if (bcf_gt_is_missing(gt[0]) || bcf_gt_is_missing(gt[1])) {
                orig_gt = GT_MISSING;

            } else { // useful

                // allow setting N/N to 1/1 later
                if (bcf_gt_allele(gt[0]) == bcf_gt_allele(gt[1])) same = true;

                if (bcf_gt_allele(gt[0]) == 0) { // REF
                    switch (bcf_gt_allele(gt[1])) {
                        case 0: orig_gt = GT_REF_REF; break;
                        case 1: orig_gt = GT_REF_ALT1; break;
                        default: orig_gt = GT_OTHER; break;
                    }
                } else if (bcf_gt_allele(gt[0]) == 1) { // ALT1
                    switch (bcf_gt_allele(gt[1])) {
                        case 0: orig_gt = GT_ALT1_REF; break;
                        case 1: orig_gt = GT_ALT1_ALT1; break;
                        case 2: orig_gt = GT_ALT1_ALT2; break;
                        default: orig_gt = GT_OTHER; break;
                    }
                } else if (bcf_gt_allele(gt[0]) == 2) { // ALT2
                    orig_gt = (bcf_gt_allele(gt[1]) == 1) ? GT_ALT2_ALT1 : GT_OTHER;
                } else {
                    orig_gt = GT_OTHER;
                }
            }

        } else if (ngt > 2) { // polyploid
            ERROR("Expected monoploid/diploid %s VCF, found variant with ploidy %d",
                    callset_strs[callset].data(), ngt);
        }
        GT_counts[orig_gt]++;

        // parse PS: https://github.com/samtools/htslib/blob/99415e2a2ce26bdbf4e910954330ea769de2c3f0/htslib/vcf.h#L1096
        nPS = bcf_get_format_int32(hdr, rec, "PS", &PS, &PS_memsize);
        if (nPS == -1) { // PS not defined in header
            if (!PS_warn) {
                PS_warn = true;
                WARN("'PS' tag not defined in %s VCF header, assuming one phase set per contig",
                        callset_strs[callset].data());
            }
            phase_set = 0;

        } else if (nPS == -3) { // PS tag missing
            // only matters if not haploid and GTs differ for this variant
            if (ngt > 1 && bcf_gt_allele(gt[0]) != bcf_gt_allele(gt[1])) {
                if (g.verbosity > 1)
                    WARN("No PS tag in %s VCF at %s:%lld",
                            callset_strs[callset].data(), ctg.data(), (long long)rec->pos);
                PS_missing_total++;
                phase_set = 0;
            }

        } else if (nPS <= 0) { // other error
            ERROR("Failed to read %s PS at %s:%lld", 
                    callset_strs[callset].data(), ctg.data(), (long long)rec->pos);
        } else {
            phase_set = PS[0];
        }


        // parse variant type
        for (int hap = 0; hap < std::abs(ngt); hap++) { // allow single-allele chrX, chrY

            // set simplified GT (0|1, 1|0, or 1|1), (0|0 and .|. skipped later)
            int simple_gt = hap ? GT_REF_ALT1 : GT_ALT1_REF; // 0|1 or 1|0 default
            if (same) simple_gt = GT_ALT1_ALT1; // overwrite 1|1 if both agree

            // get ref and allele, skipping ref query
            std::string ref = rec->d.allele[0];
            int alt_idx = ngt < 0 ? 1 : bcf_gt_allele(gt[hap]); // if no GT, assume 1
            if (alt_idx < 0) {
                if (g.verbosity > 1)
                    WARN("Variant with unknown allele (.) in %s VCF at %s:%lld, skipping",
                        callset_strs[callset].data(), ctg.data(), (long long)rec->pos);
                unknown_allele_total += 1;
                continue;
            }
            if (alt_idx == 0) continue; // nothing to do if reference
            std::string alt = rec->d.allele[alt_idx];

            // skip unphased heterozygous variants (1/1 is allowed, 0/1 is not)
            if (ngt == 2 && !same && !bcf_gt_is_phased(gt[HAP2])) { // only HAP2 is set, not sure why...
                if (g.verbosity > 1) {
                    WARN("Variant with unphased genotype in %s VCF at %s:%lld %s %s, skipping",
                        callset_strs[callset].data(), ctg.data(), (long long)rec->pos, ref.data(), alt.data());
                }
                unphased_gt_total += 1;
                continue; 
            }

            // skip spanning deletion
            if (alt == "*") { 
                ntypes[hap][TYPE_REF]++; 
                spanning_del_total++;
                continue; 
            }

            // determine variant type
            int pos = rec->pos;
            int type = -1;
            int lm = 0; // match from left->right (trim prefix)
            int rm = -1;// match from right->left (simplify complex variants CPX->INDEL)
            int reflen = int(ref.size());
            int altlen = int(alt.size());
            if (altlen-reflen > 0) { // insertion
                while (lm < reflen && ref[lm] == alt[lm]) lm++;
                while (reflen+rm >= lm && 
                        ref[reflen+rm] == alt[altlen+rm]) rm--;
                if (lm > reflen+rm) type = TYPE_INS; else type = TYPE_CPX;
                pos += lm;
                alt = alt.substr(lm, altlen+rm-lm+1);
                ref = ref.substr(lm, reflen+rm-lm+1);

            } else if (altlen-reflen < 0) { // deletion
                while (lm < altlen && ref[lm] == alt[lm]) lm++;
                while (altlen+rm >= lm && 
                        ref[reflen+rm] == alt[altlen+rm]) rm--;
                if (lm > altlen+rm) type = TYPE_DEL; else type = TYPE_CPX;
                pos += lm;
                alt = alt.substr(lm, altlen+rm-lm+1);
                ref = ref.substr(lm, reflen+rm-lm+1);

            } else { // substitution
                if (ref.size() == 1) {
                    type = (ref[0] == alt[0] ? TYPE_REF : TYPE_SUB);
                    if (type == TYPE_REF) {
                        ref_call_total++;
                        continue;
                    }
                } else {
                    if (ref.substr(1) == alt.substr(1)){
                        type = TYPE_SUB;
                        ref = ref[0]; alt = alt[0]; // chop off matches
                    }
                    else type = TYPE_CPX;
                }
            }

            // calculate reference length of variant
            int rlen = 0;
            switch (type) {
                case TYPE_INS:
                    rlen = 0; break;
                case TYPE_SUB:
                case TYPE_REF:
                    rlen = 1; break;
                case TYPE_DEL:
                case TYPE_CPX:
                    rlen = ref.size(); break;
                default:
                    ERROR("Unexpected variant type: %d", type);
                    break;
            }

            // check that variant (original representation) is in region of interest
            uint8_t loc = g.bed.contains(ctg, rec->pos, rec->pos + reflen, type);
            switch (loc) {
                case BED_OUTSIDE: 
                case BED_OFFCTG:
                case BED_BORDER:
                    nregions[loc]++;
                    continue; // discard variant
                case BED_INSIDE: 
                    nregions[loc]++;
                    break;
                default:
                    ERROR("Unexpected BED region type: %d", loc);
                    break;
            }

            // skip variants that are too large
            if (int(ref.size()) > g.max_size || int(alt.size()) > g.max_size) {
                if (g.verbosity > 1)
                    WARN("Large variant of length %d in %s VCF at %s:%lld, skipping",
                        int(std::max(ref.size(), alt.size())),
                        callset_strs[callset].data(), ctg.data(), (long long)rec->pos);
                too_large_var_total++;
                continue;
            }

            // skip overlapping variants
            if (prev_end[hap] > pos || // skip overlapping variants
                    (prev_end[hap] == pos && prev_type[hap] == TYPE_INS && type == TYPE_INS)) { // don't allow two insertions at same position
                if (g.verbosity > 1) {
                    WARN("Overlap in %s VCF variants at %s:%i, skipping", 
                            callset_strs[callset].data(), ctg.data(), pos);
                }
                overlapping_var_total++;
                continue;
            }
            // update simple_gt if corresponding variant on other hap is skipped
            if (simple_gt == GT_ALT1_ALT1 && (prev_end[hap^1] > pos ||
                    (prev_end[hap^1] == pos && prev_type[hap^1] == TYPE_INS && type == TYPE_INS))) {
                simple_gt = hap ? GT_REF_ALT1 : GT_ALT1_REF;
            }

            // add to haplotype-specific query info
            if (int(ref.size()) > g.max_cluster_var_size || int(alt.size()) > g.max_cluster_var_size) {
                large_var_total++;
                if (type == TYPE_CPX) { // split CPX into INS+DEL
                    large_variant_data->variants[hap][ctg]->add_var(pos, 0, // INS
                        TYPE_INS, loc, "", alt, simple_gt, ngq ? gq[0]:0, vq, phase_set);
                    large_variant_data->variants[hap][ctg]->add_var(pos, rlen, // DEL
                        TYPE_DEL, loc, ref, "", simple_gt, ngq ? gq[0]:0, vq, phase_set);
                    complex_total++;
                } else {
                    large_variant_data->variants[hap][ctg]->add_var(pos, rlen,
                            type, loc, ref, alt, simple_gt, ngq ? gq[0]:0, vq, phase_set);
                }
            } else {
                if (type == TYPE_CPX) { // split CPX into INS+DEL
                    variant_data->variants[hap][ctg]->add_var(pos, 0, // INS
                        TYPE_INS, loc, "", alt, simple_gt, ngq ? gq[0]:0, vq, phase_set);
                    variant_data->variants[hap][ctg]->add_var(pos, rlen, // DEL
                        TYPE_DEL, loc, ref, "", simple_gt, ngq ? gq[0]:0, vq, phase_set);
                    complex_total++;
                } else {
                    variant_data->variants[hap][ctg]->add_var(pos, rlen,
                            type, loc, ref, alt, simple_gt, ngq ? gq[0]:0, vq, phase_set);
                }
            }

            prev_end[hap] = pos + rlen;
            prev_type[hap] = type;
            npass[hap]++;
            ntypes[hap][type]++;
        }
    }

    // SUMMARY PRINTING

    /* if (gq_missing_total) */ 
    /*     WARN("%d total missing GQ tags in %s VCF, all considered GQ=0", */
    /*         gq_missing_total, callset_strs[callset].data()); */

    if (failed_filter_total && print) 
        INFO("%d variants failed FILTER in %s VCF, skipped",
            failed_filter_total, callset_strs[callset].data());

    if (pass_min_qual[FALSE] && print)
        INFO("%d variants of low quality (<%d) in %s VCF, skipped", 
            pass_min_qual[FALSE], g.min_qual, callset_strs[callset].data());

    if (wrong_ploidy_total) 
        WARN("%d variants with incorrect ploidy in %s VCF, kept",
            wrong_ploidy_total, callset_strs[callset].data());

    if (print) INFO("  Genotypes:");
    for (size_t i = 0; i < gt_strs.size(); i++) {
        if (print && GT_counts[i]) INFO("    %3s: %i", gt_strs[i].data(), GT_counts[i]);
    }
    if (float(GT_counts[GT_REF_ALT1]) / (GT_counts[GT_ALT1_REF]+1) > 2 ||
        float(GT_counts[GT_ALT1_REF]) / (GT_counts[GT_REF_ALT1]+1) > 2)
        WARN("Imbalance of heterozygous variant phasing, VCF may be improperly phased")
    if (print) INFO(" ");

    if (PS_missing_total) 
        WARN("%d variants missing PS tags in %s VCF, kept",
            PS_missing_total, callset_strs[callset].data());

    multi_total = GT_counts[GT_ALT1_ALT1] + GT_counts[GT_ALT1_ALT2] +
        GT_counts[GT_ALT2_ALT1] + GT_counts[GT_OTHER];
    if (multi_total && print)
        INFO("%d homozygous and multi-allelic variants in %s VCF, split for evaluation",
            multi_total, callset_strs[callset].data());

    if (unknown_allele_total) 
        WARN("%d variants with unknown (./.) alleles in %s VCF, skipped",
            unknown_allele_total, callset_strs[callset].data());

    if (unphased_gt_total) 
        WARN("%d variants with unphased genotypes in %s VCF, skipped",
            unphased_gt_total, callset_strs[callset].data());

    if (spanning_del_total)
        WARN("%d variants spanned by deletion in %s VCF, skipped", 
                spanning_del_total, callset_strs[callset].data());

    if (ref_call_total) 
        WARN("%d reference variants in %s VCF, skipped",
            ref_call_total, callset_strs[callset].data());

    if (nregions[BED_OFFCTG] + nregions[BED_OUTSIDE] && print)
        INFO("%d variants outside selected regions in %s VCF, skipped",
                nregions[BED_OFFCTG] + nregions[BED_OUTSIDE], 
                callset_strs[callset].data());

    if (nregions[BED_BORDER] && print)
        INFO("%d variants on border of selected regions in %s VCF, skipped",
                nregions[BED_BORDER], callset_strs[callset].data());

    if (too_large_var_total)
        WARN("%d large (size > %d) variants in %s VCF, skipped", 
                too_large_var_total, g.max_size, callset_strs[callset].data());

    if (large_var_total)
        INFO("%d large (size > %d) variants in %s VCF, set aside for second round evaluation", 
                large_var_total, g.max_cluster_var_size, callset_strs[callset].data());

    if (overlapping_var_total)
        WARN("%d overlapping variants in %s VCF, skipped", 
                overlapping_var_total, callset_strs[callset].data());

    if (complex_total)
        WARN("%d complex (CPX) variants in %s VCF, split into INS + DEL", 
                complex_total, callset_strs[callset].data());

    if (print) INFO(" ");
    if (print) INFO("  Variant types:");
    if (g.verbosity >= 2) { // show each hap separately
        for (int h = 0; h < HAPS; h++) {
            if (print) INFO("    Haplotype %i", h+1);
            for (size_t i = 0; i < type_strs.size(); i++) {
                if (print) INFO("      %s: %i", type_strs[i].data(), ntypes[h][i]);
            }
        }
        if (print) INFO(" ");
    } else { // summarize
        for (size_t i = 0; i < type_strs.size(); i++) {
            if (print && ntypes[HAP1][i] + ntypes[HAP2][i]) 
                INFO("    %s: %i", type_strs[i].data(), ntypes[HAP1][i] + ntypes[HAP2][i]);
        }
    }
    if (print) INFO(" ");

    if (print) INFO("  Contigs:");
    for (size_t i = 0; i < variant_data->contigs.size(); i++) {
        if (print) INFO("    [%2lu] %s: %d | %d variants", i, variant_data->contigs[i].data(),
                variant_data->variants[HAP1][variant_data->contigs[i]]->n, 
                variant_data->variants[HAP2][variant_data->contigs[i]]->n);
    }
    if (print) INFO(" ");

    if (print) INFO("  %s VCF overview:", callset_strs[callset].data());
    if (print) INFO("      Total %s variants: %d", 
            callset_strs[callset].data(),  n + multi_total + complex_total);
    if (print) INFO("      Kept  %s variants: %d", 
            callset_strs[callset].data(), npass[HAP1] + npass[HAP2] + complex_total);

    free(gq);
    free(fgq);
    free(gt);
    free(PS);
    free(ctgnames);
    bcf_hdr_destroy(hdr);
    bcf_close(vcf);
    bcf_destroy(rec);
    return;
error2:
    free(ctgnames);
error1:
    bcf_close(vcf);
    bcf_hdr_destroy(hdr);
    return;

}
