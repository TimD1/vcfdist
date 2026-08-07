/**
 * @file variant.cpp
 * @brief Per-contig and per-callset variant containers with VCF parsing and output utilities.
 */
#include <algorithm>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <vector>
#include <cmath>

#include "htslib/vcf.h"

#include "variant.h"
#include "print.h"
#include "dist.h"

/**************************************************************************************************/

/**
 * @brief Constructs a contig-specific variant container with empty data structures.
 * @param[in] ctg Contig name
 */
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


/**************************************************************************************************/


/**
 * @brief Appends a variant with all fields explicitly specified.
 * @param[in] var Every field describing the variant; gt_qual and var_qual are capped at g.max_qual
 */
// TODO: remove assumption that no variants match
void ctgVariants::add_var(const var_fields & var) {
    // set for all variants
    this->poss.push_back(var.pos);
    this->rlens.push_back(var.rlen);
    this->types.push_back(var.type);
    this->locs.push_back(var.loc);
    this->refs.push_back(var.ref);
    this->alts.push_back(var.alt);
    this->orig_gts.push_back(var.orig_gt);
    this->gt_quals.push_back(std::min(var.gt_qual, float(g.max_qual)));
    this->var_quals.push_back(std::min(var.var_qual, float(g.max_qual)));
    this->phase_sets.push_back(var.phase_set);
    this->rec_idxs.push_back(var.rec_idx);
    this->alt_idxs.push_back(var.alt_idx);
    this->ploidies.push_back(var.ploidy);
    this->superclusters.push_back(var.supercluster);
    this->n++;

    // added during precision/recall analysis
    this->calc_gts.push_back(var.calc_gt);
    for (int hap = 0; hap < HAPS; hap++) {
        this->errtypes[hap].push_back(var.hap[hap].errtype);
        this->sync_group[hap].push_back(var.hap[hap].sync_group);
        this->callq[hap].push_back(var.hap[hap].callq);
        this->ref_ed[hap].push_back(var.hap[hap].ref_ed);
        this->query_ed[hap].push_back(var.hap[hap].query_ed);
        this->credit[hap].push_back(var.hap[hap].credit);
    }

    // added during phasing analysis
    this->phases.push_back(PHASE_NONE);
    this->pb_phases.push_back(PHASE_NONE);
    this->ac_errtype.push_back(AC_UNKNOWN);
}


/**************************************************************************************************/


/**
 * @brief Returns every field of one variant, for copying it into another container.
 * @param[in] idx Variant index
 * @return All fields of variant idx, suitable for passing to another container's add_var()
 * @throws ERROR Variant index out of range for this container
 */
var_fields ctgVariants::get_var(int idx) const {
    if (idx < 0 || idx >= this->n) {
        ERROR("Variant index %d out of range for contig %s (%d variants) in get_var()",
                idx, this->ctg.data(), this->n);
    }
    var_fields var = {
        .pos = this->poss[idx],
        .rlen = this->rlens[idx],
        .type = this->types[idx],
        .loc = this->locs[idx],
        .ref = this->refs[idx],
        .alt = this->alts[idx],
        .orig_gt = this->orig_gts[idx],
        .gt_qual = this->gt_quals[idx],
        .var_qual = this->var_quals[idx],
        .phase_set = this->phase_sets[idx],
        .rec_idx = this->rec_idxs[idx],
        .alt_idx = this->alt_idxs[idx],
        .ploidy = this->ploidies[idx],
        .supercluster = this->superclusters[idx],
        .calc_gt = this->calc_gts[idx],
    };
    for (int hap = 0; hap < HAPS; hap++) {
        var.hap[hap] = {
            .errtype = this->errtypes[hap][idx],
            .sync_group = this->sync_group[hap][idx],
            .callq = this->callq[hap][idx],
            .ref_ed = this->ref_ed[hap][idx],
            .query_ed = this->query_ed[hap][idx],
            .credit = this->credit[hap][idx],
        };
    }
    return var;
}


/**************************************************************************************************/


/**
 * @brief Classifies variant as SNP, INDEL, or SV based on reference length and g.sv_threshold.
 * @param[in] vi Variant index
 * @return VARTYPE_SNP, VARTYPE_INDEL, or VARTYPE_SV
 */
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


/**************************************************************************************************/


/**
 * @brief Returns the alternate allele count of a diploid genotype.
 * @param[in] gt Genotype (GT_*)
 * @return 0, 1, or 2 alternate alleles, or -1 if gt carries no diploid allele count
 */
static int allele_count(uint8_t gt) {
    switch (gt) {
        case GT_REF_REF:   return 0;
        case GT_REF_ALT1:
        case GT_ALT1_REF:  return 1;
        case GT_ALT1_ALT1: return 2;
        default:           return -1;
    }
}


/**
 * @brief Maps a site's truth and query alternate allele counts to an allele count error type.
 * @param[in] truth_ac Truth alternate allele count, or -1 if unknown
 * @param[in] query_ac Query alternate allele count, or -1 if unknown
 * @return AC_ERR_*_TO_*, or AC_UNKNOWN if either count is unknown or both are zero
 */
static int ac_errtype_from_counts(int truth_ac, int query_ac) {
    switch (truth_ac) {
        case 0:
            if (query_ac == 1) return AC_ERR_0_TO_1;
            if (query_ac == 2) return AC_ERR_0_TO_2;
            break;
        case 1:
            if (query_ac == 0) return AC_ERR_1_TO_0;
            if (query_ac == 1) return AC_ERR_1_TO_1;
            if (query_ac == 2) return AC_ERR_1_TO_2;
            break;
        case 2:
            if (query_ac == 0) return AC_ERR_2_TO_0;
            if (query_ac == 1) return AC_ERR_2_TO_1;
            if (query_ac == 2) return AC_ERR_2_TO_2;
            break;
    }
    return AC_UNKNOWN;
}


/**
 * @brief Records a variant's allele count error type from its original and calculated genotypes.
 *
 * The value keeps one absolute truth-allele-count-then-query-allele-count direction on both
 * callsets, so the two records of a matched site report it identically. Which genotype supplies
 * which count is what differs: a record's orig_gt is its own callset's call and its calc_gt is the
 * other callset's genotype as recovered by alignment.
 *
 * @param[in] vi Variant index
 * @param[in] query True if this container holds query variants, false for truth variants
 * @return Allele count error type (AC_ERR_*_TO_* or AC_UNKNOWN); also stored in ac_errtype[vi]
 */
int ctgVariants::set_allele_errtype(int vi, bool query) {
    int truth_ac = allele_count(query ? this->calc_gts[vi] : this->orig_gts[vi]);
    int query_ac = allele_count(query ? this->orig_gts[vi] : this->calc_gts[vi]);
    return this->ac_errtype[vi] = ac_errtype_from_counts(truth_ac, query_ac);
}


/**************************************************************************************************/


/**
 * @brief Returns true if haplotypes should be swapped when reporting calc_gt data relative to orig_gt.
 * @param[in] vi Variant index
 * @return False for matching genotypes, homozygous calls, or when calc_gt is 0/0
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

/**
 * @brief Returns true if a variant is present on the specified haplotype.
 * @param[in] var_idx Variant index
 * @param[in] hap Haplotype index (0 or 1)
 * @param[in] calc If true, check calc_gts; if false, check orig_gts
 * @return True if variant is on the specified haplotype
 */
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

/**************************************************************************************************/

/**
 * @brief Sets or unsets the alternate allele on one haplotype for a calculated genotype.
 * @param[in] var_idx Variant index
 * @param[in] hap Haplotype index (0 or 1)
 * @param[in] set If true, set alternate; if false, unset it
 * @param[in] ignore_errors If true, suppress error messages for invalid transitions
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

/**************************************************************************************************/

/**
 * @brief Writes fixed VCF fields (CHROM, POS, ID, REF, ALT, QUAL, FILTER, INFO, FORMAT) for one variant.
 * @param[in] out_fp Open file pointer to output VCF
 * @param[in] ref Reference FASTA data for retrieving flanking bases for indels
 * @param[in] ctg Contig name
 * @param[in] idx Variant index in this container
 * @throws ERROR An INS/DEL sits at the contig start (0-based pos 0), leaving no preceding base to anchor
 * @throws ERROR The variant type is not TYPE_SUB, TYPE_INS, or TYPE_DEL
 */
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
        // INS/DEL are left-anchored on the preceding reference base; at contig start (0-based
        // pos 0) there is no preceding base, so guard against the out-of-bounds read of index -1
        if (this->poss[idx] == 0)
            ERROR("Cannot left-anchor INS/DEL at contig start (0-based pos 0) on '%s' in print_var_info",
                    ctg.data());
        ref_base = ref->fasta.at(ctg)[this->poss[idx]-1];
        fprintf(out_fp, "%s\t%d\t.\t%s\t%s\t.\tPASS\t.\tGT:BD:BC:RD:QD:BK:QQ:SC:SG:PS:PB:BS:VP:FE:GE", ctg.data(), 
                this->poss[idx], (ref_base + this->refs[idx]).data(), 
                (ref_base + this->alts[idx]).data());
        break;
    default:
        ERROR("print_var_info not implemented for type %d", this->types[idx]);
    }
}


/**
 * @brief Writes dot-separated empty sample fields for a variant with no call on this haplotype.
 * @param[in] out_fp Open file pointer to output VCF
 * @param[in] sc_idx Supercluster index for SC field
 * @param[in] phase_block Phase block index for PB field
 * @param[in] query If true, append newline (end of record); if false, tab (more samples follow)
 */
void ctgVariants::print_var_empty(FILE* out_fp, int sc_idx,
        int phase_block, bool query /* = false */) {
    fprintf(out_fp, "\t.:.:.:.:.:.:.:%d:.:.:%d:.:.:.:.%s", sc_idx, phase_block, query ? "\n" : "");
}


/**
 * @brief Writes sample-specific FORMAT fields for one variant to output VCF.
 * @param[in] out_fp Open file pointer to output VCF
 * @param[in] vi Variant index in this container
 * @param[in] hi Haplotype index (0 or 1)
 * @param[in] gt Genotype string (e.g., "0|1", "1|1")
 * @param[in] sc_idx Supercluster index for SC field
 * @param[in] phase_block Phase block index for PB field
 * @param[in] phase_switch True if phase switched at this position
 * @param[in] phase_flip True if phase flipped (error) at this position
 * @param[in] query If true, format as query sample; if false, as truth sample
 */
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


/**************************************************************************************************/


/**
 * @brief Parses a CIGAR string and adds resulting variants to the container.
 * @param[in] cigar CIGAR operation vector (alternating operation codes and lengths)
 * @param[in] hap Haplotype index (0 or 1)
 * @param[in] ref_pos Starting reference position
 * @param[in] ctg Contig name
 * @param[in] query Query sequence
 * @param[in] ref Reference sequence
 * @param[in] qual Quality score for all variants
 * @param[in] phase_set Phase set identifier for all variants
 * @throws ERROR Unexpected CIGAR/pointer operation not in {PTR_MAT, PTR_SUB, PTR_DEL, PTR_INS}
 */
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
                this->variants[hap][ctg]->add_var(var_fields{.pos = ref_pos+ref_idx, .rlen = 1,
                        .type = TYPE_SUB, .loc = BED_INSIDE,
                        .ref = std::string(1,ref[ref_idx]),
                        .alt = std::string(1,query[query_idx]),
                        .orig_gt = GT_REF_REF, .gt_qual = float(g.max_qual),
                        .var_qual = float(qual), .phase_set = phase_set});
                ref_idx++;
                query_idx++;
                break;

            case PTR_DEL: // deletion
                cig_idx++; indel_len++;

                // multi-base deletion
                while (cig_idx < cigar.size() && cigar[cig_idx] == PTR_DEL) {
                    cig_idx++; indel_len++;
                }
                this->variants[hap][ctg]->add_var(var_fields{.pos = ref_pos+ref_idx,
                        .rlen = indel_len, .type = TYPE_DEL, .loc = BED_INSIDE,
                        .ref = ref.substr(ref_idx, indel_len), .alt = "",
                        .orig_gt = GT_REF_REF, .gt_qual = float(g.max_qual),
                        .var_qual = float(qual), .phase_set = phase_set});
                ref_idx += indel_len;
                break;

            case PTR_INS: // insertion
                cig_idx++; indel_len++;

                // multi-base insertion
                while (cig_idx < cigar.size() && cigar[cig_idx] == PTR_INS) {
                    cig_idx++; indel_len++;
                }
                this->variants[hap][ctg]->add_var(var_fields{.pos = ref_pos+ref_idx,
                        .rlen = 0, .type = TYPE_INS, .loc = BED_INSIDE, .ref = "",
                        .alt = query.substr(query_idx, indel_len),
                        .orig_gt = GT_REF_REF, .gt_qual = float(g.max_qual),
                        .var_qual = float(qual), .phase_set = phase_set});
                query_idx += indel_len;
                break;

            default:
                ERROR("Unexpected CIGAR operation (%d) in add_variants", cigar[cig_idx]);
        }
    }
}

/**************************************************************************************************/

/**
 * @brief Constructs an empty variant data container defaulting to QUERY callset.
 */
variantData::variantData() : callset(QUERY), variants(HAPS) { ; }

/**
 * @brief Parses variants from a VCF file into a variantData container, with filtering and validation.
 * @param[in] vcf_fn Input VCF filename
 * @param[out] variant_data Container to populate with parsed variants
 * @param[in] reference Reference FASTA data for coordinate validation
 * @param[in] callset QUERY or TRUTH callset identifier
 * @throws Various errors for malformed VCF or invalid reference coordinates
 * @throws WARNING Per-reason summary totals for records dropped or altered at parse time: no-call
 *         and half-call genotypes, unphased heterozygous genotypes, spanning deletions, reference
 *         calls, missing PS tags, oversized variants, overlapping variants, and complex variants
 *         split into INS + DEL
 */
void parse_variants(const std::string & vcf_fn,
        std::shared_ptr<variantData> variant_data,
        std::shared_ptr<fastaData> reference,
        int callset) {

    // set reference fasta pointer
    variant_data->ref = reference;
    variant_data->filename = vcf_fn;

    if (callset < 0 || callset >= CALLSETS)
        ERROR("Invalid callset (%d).", callset);
    variant_data->callset = callset;

    if (g.verbosity >= 1) INFO(" ");
    if (g.verbosity >= 1) INFO("%s[%s %d/%d] Parsing %s VCF%s '%s'", COLOR_PURPLE,
            callset == QUERY ? "Q" : "T", int(idx(TIME_READ)), int(idx(TIME_TOTAL))-1,
            callset_strs[callset].data(), 
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
    std::unordered_set<int> prev_rids; // contigs already parsed, to reject an unsorted VCF
    int prev_rid = -1;
    std::unordered_map<int, int> ctglens;
    std::string ctg;
    EnumArray<bedloc_t, int, BEDLOC_SLOTS> nregions{};
    std::vector<int> pass_min_qual(2, 0);

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

    /* int gq_missing_total = 0; */
    int PS_missing_total = 0;
    int overlapping_var_total = 0;
    int spanning_del_total = 0;
    int unknown_allele_total = 0; // no-call records (.|. or .), dropped entirely
    int half_call_total = 0;      // half-call records (1|. or .|1), known allele kept
    int unphased_gt_total = 0;
    int too_large_var_total = 0;
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
                prev_rids.insert(rec->rid);
                variant_data->contigs.push_back(ctg);
                variant_data->observed_ploidies.push_back({});
                variant_data->lengths.push_back(ctglens[rec->rid]);
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

        // record this record's ploidy; mixed ploidy within a contig is legitimate, as on a chrX
        // carrying both PAR (diploid) and non-PAR (haploid) calls, so nothing is enforced here
        int ctg_idx = std::find(variant_data->contigs.begin(), variant_data->contigs.end(), ctg)
                - variant_data->contigs.begin();
        variant_data->observed_ploidies[ctg_idx].insert(std::abs(ngt));

        // parse genotype info
        int orig_gt = GT_REF_REF;
        bool same = false;
        if (ngt == -1) { // no info, assume monoploid
            orig_gt = GT_ALT1;
        } else if (ngt == 1) { // monoploid/haploid

            // set 1 if allele_idx > 0, no-call if the single allele is missing
            if (bcf_gt_is_missing(gt[0])) orig_gt = GT_MISSING;
            else orig_gt = bcf_gt_allele(gt[0]) ? GT_ALT1 : GT_REF;

        } else if (ngt == 2) { // diploid

            // distinguish a no-call (both alleles missing) from a half call (exactly one missing)
            bool hap1_missing = bcf_gt_is_missing(gt[HAP1]);
            bool hap2_missing = bcf_gt_is_missing(gt[HAP2]);
            if (hap1_missing && hap2_missing) { // no call (.|.), record is dropped
                orig_gt = GT_MISSING;

            } else if (hap1_missing || hap2_missing) { // half call (1|.), known allele is kept
                orig_gt = GT_HALF;

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

        // count missing alleles once per record, not once per haplotype
        if (orig_gt == GT_MISSING) {
            if (g.verbosity > 1)
                WARN("Variant with no known alleles (.|.) in %s VCF at %s:%lld, skipping",
                    callset_strs[callset].data(), ctg.data(), (long long)rec->pos);
            unknown_allele_total += 1;
        } else if (orig_gt == GT_HALF) {
            if (g.verbosity > 1)
                WARN("Variant with a half call (1|.) in %s VCF at %s:%lld, keeping known allele",
                    callset_strs[callset].data(), ctg.data(), (long long)rec->pos);
            half_call_total += 1;
        }

        // parse PS: https://github.com/samtools/htslib/blob/99415e2a2ce26bdbf4e910954330ea769de2c3f0/htslib/vcf.h#L1096
        // scoped to this record, so one lacking a PS tag can never inherit the previous record's
        int phase_set = 0; // no declared phase set: one implicit phase set per contig
        nPS = bcf_get_format_int32(hdr, rec, "PS", &PS, &PS_memsize);
        if (nPS == -1) { // PS not defined in header
            if (!PS_warn) {
                PS_warn = true;
                WARN("'PS' tag not defined in %s VCF header, assuming one phase set per contig",
                        callset_strs[callset].data());
            }

        } else if (nPS == -3) { // PS tag missing
            // only counted if not haploid and GTs differ, since only then is phase unresolved
            if (ngt > 1 && bcf_gt_allele(gt[0]) != bcf_gt_allele(gt[1])) {
                if (g.verbosity > 1)
                    WARN("No PS tag in %s VCF at %s:%lld",
                            callset_strs[callset].data(), ctg.data(), (long long)rec->pos);
                PS_missing_total++;
            }

        } else if (nPS <= 0) { // other error
            ERROR("Failed to read %s PS at %s:%lld", 
                    callset_strs[callset].data(), ctg.data(), (long long)rec->pos);
        } else {
            phase_set = PS[0];
        }


        // snapshot each hap's previous variant, before this record's own copies overwrite it below
        std::vector<int> rec_prev_end = prev_end;
        std::vector<int> rec_prev_type = prev_type;

        // parse variant type
        for (int hap = 0; hap < std::abs(ngt); hap++) { // allow single-allele chrX, chrY

            // set simplified GT (0|1, 1|0, or 1|1), (0|0 and .|. skipped later)
            int simple_gt = hap ? GT_REF_ALT1 : GT_ALT1_REF; // 0|1 or 1|0 default
            if (same) simple_gt = GT_ALT1_ALT1; // overwrite 1|1 if both agree

            // get ref and allele, skipping ref query
            std::string ref = rec->d.allele[0];
            int alt_idx = ngt < 0 ? 1 : bcf_gt_allele(gt[hap]); // if no GT, assume 1
            if (alt_idx < 0) continue; // missing allele (.), counted once per record above
            if (alt_idx == 0) continue; // nothing to do if reference
            std::string alt = rec->d.allele[alt_idx];

            // uppercase before any comparison, since soft-masked reference sequence reaches us as
            // lowercase and case must not decide whether alleles match
            std::transform(ref.begin(), ref.end(), ref.begin(), ::toupper);
            std::transform(alt.begin(), alt.end(), alt.begin(), ::toupper);

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

                // skip reference calls, where ALT is identical to REF (e.g. A -> A, AT -> AT)
                if (ref == alt) {
                    ref_call_total++;
                    continue;
                }
                if (ref.size() == 1) {
                    type = TYPE_SUB;
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
            bedloc_t loc = g.bed.contains(ctg, rec->pos, rec->pos + reflen, type);
            switch (loc) {
                case BED_OUTSIDE: 
                case BED_OFFCTG:
                case BED_BORDER:
                    nregions[loc]++;
                    continue; // discard variant
                case BED_INSIDE:
                    nregions[loc]++;
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
            if (simple_gt == GT_ALT1_ALT1 && (rec_prev_end[hap^1] > pos ||
                    (rec_prev_end[hap^1] == pos && rec_prev_type[hap^1] == TYPE_INS &&
                     type == TYPE_INS))) {
                simple_gt = hap ? GT_REF_ALT1 : GT_ALT1_REF;
            }

            // add to haplotype-specific query info
            int rec_idx = n - 1; // 0-based ordinal of this record within the input VCF
            uint8_t ploidy = uint8_t(std::abs(ngt));
            // both CPX halves derive from the same original allele, so they share alt_idx
            if (type == TYPE_CPX) { // split CPX into INS+DEL
                variant_data->variants[hap][ctg]->add_var(var_fields{.pos = pos, .rlen = 0, // INS
                    .type = TYPE_INS, .loc = loc, .ref = "", .alt = alt,
                    .orig_gt = uint8_t(simple_gt), .gt_qual = float(ngq ? gq[0]:0),
                    .var_qual = vq, .phase_set = phase_set,
                    .rec_idx = rec_idx, .alt_idx = alt_idx, .ploidy = ploidy});
                variant_data->variants[hap][ctg]->add_var(var_fields{.pos = pos, .rlen = rlen, // DEL
                    .type = TYPE_DEL, .loc = loc, .ref = ref, .alt = "",
                    .orig_gt = uint8_t(simple_gt), .gt_qual = float(ngq ? gq[0]:0),
                    .var_qual = vq, .phase_set = phase_set,
                    .rec_idx = rec_idx, .alt_idx = alt_idx, .ploidy = ploidy});
                complex_total++;
            } else {
                variant_data->variants[hap][ctg]->add_var(var_fields{.pos = pos, .rlen = rlen,
                        .type = uint8_t(type), .loc = loc, .ref = ref, .alt = alt,
                        .orig_gt = uint8_t(simple_gt), .gt_qual = float(ngq ? gq[0]:0),
                        .var_qual = vq, .phase_set = phase_set,
                        .rec_idx = rec_idx, .alt_idx = alt_idx, .ploidy = ploidy});
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

    if (pass_min_qual[false] && print)
        INFO("%d variants of low quality (<%d) in %s VCF, skipped",
            pass_min_qual[false], g.min_qual, callset_strs[callset].data());

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
        WARN("%d variants with no known alleles (.|.) in %s VCF, skipped",
            unknown_allele_total, callset_strs[callset].data());

    if (half_call_total)
        WARN("%d variants with a half call (1|.) in %s VCF, known allele kept",
            half_call_total, callset_strs[callset].data());

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
