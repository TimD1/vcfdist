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
#include <cstdio>

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
    this->matched_gts.push_back(var.matched_gt);
    for (hap_t hap : EnumRange<hap_t, HAP_SLOTS>{}) {
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
        .matched_gt = this->matched_gts[idx],
    };
    for (hap_t hap : EnumRange<hap_t, HAP_SLOTS>{}) {
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
sizeclass_t ctgVariants::get_vartype(int vi) {
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
 * @param[in] gt Genotype
 * @return 0, 1, or 2 alternate alleles
 * @throws ERROR A gt_t outside the four evaluation genotypes, which the type cannot represent
 */
static int allele_count(gt_t gt) {
    switch (gt) {
        case GT_REF_REF: return 0;
        case GT_REF_ALT:
        case GT_ALT_REF: return 1;
        case GT_ALT_ALT: return 2;
    }
    ERROR("Unexpected genotype %d in allele_count()", int(gt));
}


/**
 * @brief Maps a site's truth and query alternate allele counts to an allele count error type.
 * @param[in] truth_ac Truth alternate allele count (0, 1, or 2)
 * @param[in] query_ac Query alternate allele count (0, 1, or 2)
 * @return AC_ERR_*_TO_*, or AC_UNKNOWN if both counts are zero
 */
static ac_errtype_t ac_errtype_from_counts(int truth_ac, int query_ac) {
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
 * @brief Records a variant's allele count error type from its original and matched genotypes.
 *
 * The value keeps one absolute truth-allele-count-then-query-allele-count direction on both
 * callsets, so the two records of a matched site report it identically. Which genotype supplies
 * which count is what differs: a record's orig_gt is its own callset's call and its matched_gt is the
 * other callset's genotype as recovered by alignment.
 *
 * @param[in] vi Variant index
 * @param[in] query True if this container holds query variants, false for truth variants
 * @return Allele count error type (AC_ERR_*_TO_* or AC_UNKNOWN); also stored in ac_errtype[vi]
 */
ac_errtype_t ctgVariants::set_allele_errtype(int vi, bool query) {
    int truth_ac = allele_count(query ? this->matched_gts[vi] : this->orig_gts[vi]);
    int query_ac = allele_count(query ? this->orig_gts[vi] : this->matched_gts[vi]);
    return this->ac_errtype[vi] = ac_errtype_from_counts(truth_ac, query_ac);
}


/**************************************************************************************************/


/**
 * @brief Returns true if haplotypes should be swapped when reporting matched_gt data relative to orig_gt.
 * @param[in] vi Variant index
 * @return False for matching genotypes, homozygous calls, or when matched_gt is 0/0
 */
bool ctgVariants::matched_gt_is_swapped(int vi /* variant index */) const {
    // 0|0,0|0 and 0|1,0|1 and 1|0,1|0 and 1|1,1|1
    if (this->orig_gts[vi] == this->matched_gts[vi]) {
        return false;
    }
    // orig_gt == 1|1 or 0|0 or matched_gt == 0|0, all/no matched_gt data will be reported, order doesn't matter
    else if (this->orig_gts[vi] == GT_ALT_ALT ||
             this->orig_gts[vi] == GT_REF_REF ||
             this->matched_gts[vi] == GT_REF_REF) {
        return false;
    }
    // 0|1,1|0 and 1|0,0|1
    else if ((this->orig_gts[vi] == GT_REF_ALT && this->matched_gts[vi] == GT_ALT_REF) ||
             (this->orig_gts[vi] == GT_ALT_REF && this->matched_gts[vi] == GT_REF_ALT)) {
        return true;
    }
    // orig_gt = 0|1, choose better matched_gt
    else if (this->orig_gts[vi] == GT_REF_ALT && this->matched_gts[vi] == GT_ALT_ALT) {
        return this->credit[HAP1][vi] > this->credit[HAP2][vi];
    }
    // orig_gt = 1|0, choose better matched_gt
    else if (this->orig_gts[vi] == GT_ALT_REF && this->matched_gts[vi] == GT_ALT_ALT) {
        return this->credit[HAP2][vi] > this->credit[HAP1][vi];
    // the branches above are exhaustive over the four evaluation genotypes, so this is a guard
    } else {
        ERROR("Unexpected orig/matched genotypes for variant (%s -> %s) at pos %d: orig=%s matched=%s",
                this->refs[vi].data(),
                this->alts[vi].data(),
                this->poss[vi],
                gt_strs[this->orig_gts[vi]].data(),
                gt_strs[this->matched_gts[vi]].data()
        );
    }
}

/**
 * @brief Returns true if a variant is present on the specified haplotype.
 * @param[in] var_idx Variant index
 * @param[in] hap Haplotype index (0 or 1)
 * @param[in] matched If true, check matched_gts; if false, check orig_gts
 * @return True if variant is on the specified haplotype
 */
bool ctgVariants::var_on_hap(int var_idx, hap_t hap, bool matched) const {
    // simple gt, always (0|1, 1|0, or 1|1)
    gt_t gt = matched ? this->matched_gts[var_idx] : this->orig_gts[var_idx];
    if (hap == HAP1 && (gt == GT_ALT_REF || gt == GT_ALT_ALT)) return true;
    if (hap == HAP2 && (gt == GT_REF_ALT || gt == GT_ALT_ALT)) return true;
    return false;
}

/**************************************************************************************************/

/**
 * @brief Sets or unsets the alternate allele on one haplotype for a matched genotype.
 * @param[in] var_idx Variant index
 * @param[in] hap Haplotype index (0 or 1)
 * @param[in] set If true, set alternate; if false, unset it
 * @param[in] ignore_errors If true, suppress error messages for invalid transitions
 */
void ctgVariants::set_var_matched_gt_on_hap(int var_idx, hap_t hap, bool set,
        bool ignore_errors) {
    if (this->matched_gts[var_idx] == GT_REF_REF) {
        if (set) {
            this->matched_gts[var_idx] = hap == HAP1 ? GT_ALT_REF : GT_REF_ALT;
        } else { // unset
            if (!ignore_errors) ERROR("Variant matched_gt already unset for variant %d at %s:%d hap %d",
                    var_idx, this->ctg.data(), this->poss[var_idx], int(idx(hap)));
        }

    } else if (this->matched_gts[var_idx] == GT_REF_ALT) {
        if (set) {
            if (hap == HAP2) {
                if (!ignore_errors) ERROR("Variant matched_gt already set for variant %d at %s:%d hap %d",
                    var_idx, this->ctg.data(), this->poss[var_idx], int(idx(hap)));
            } else {
                this->matched_gts[var_idx] = GT_ALT_ALT;
            }
        } else { // unset
            if (hap == HAP1) {
                if (!ignore_errors) ERROR("Variant matched_gt already unset for variant %d at %s:%d hap %d",
                    var_idx, this->ctg.data(), this->poss[var_idx], int(idx(hap)));
            } else {
                this->matched_gts[var_idx] = GT_REF_REF;
            }
        }

    } else if (this->matched_gts[var_idx] == GT_ALT_REF) {
        if (set) {
            if (hap == HAP1) {
                if (!ignore_errors) ERROR("Variant matched_gt already set for variant %d at %s:%d hap %d",
                    var_idx, this->ctg.data(), this->poss[var_idx], int(idx(hap)));
            } else {
                this->matched_gts[var_idx] = GT_ALT_ALT;
            }
        } else { // unset
            if (hap == HAP2) {
                if (!ignore_errors) ERROR("Variant matched_gt already unset for variant %d at %s:%d hap %d",
                    var_idx, this->ctg.data(), this->poss[var_idx], int(idx(hap)));
            } else {
                this->matched_gts[var_idx] = GT_REF_REF;
            }
        }

    } else if (this->matched_gts[var_idx] == GT_ALT_ALT) {
        if (set) {
            if (!ignore_errors) ERROR("Variant matched_gt already set for variant %d at %s:%d hap %d",
                    var_idx, this->ctg.data(), this->poss[var_idx], int(idx(hap)));
        } else {
            this->matched_gts[var_idx] = hap == HAP1 ? GT_REF_ALT : GT_ALT_REF;
        }

    // the branches above are exhaustive over the four evaluation genotypes, so this is a guard
    } else {
        ERROR("Unexpected matched_gts value '%s' in set_var_matched_gt_on_hap() for variant %d at %s:%d",
            gt_strs[this->matched_gts[var_idx]].data(), var_idx, this->ctg.data(), this->poss[var_idx]);
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
        ERROR("print_var_info not implemented for type %d", static_cast<int>(this->types[idx]));
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
 * @brief Renders a credit exactly as printf's "%f" would, for embedding in a comma-separated list.
 * @param[in] credit Credit on the interval [0,1]
 * @return The credit with six digits after the decimal point
 */
static std::string credit_str(float credit) {
    char buf[32];
    snprintf(buf, sizeof(buf), "%f", credit);
    return std::string(buf);
}


/**
 * @brief Writes sample-specific FORMAT fields for one variant to output VCF.
 *
 * One record is written per variant rather than per haplotype, so the per-haplotype fields (BD, BC,
 * RD, QD, BK, SG) are comma-separated lists carrying one value per allele of the emitted GT. GT is
 * rendered from orig_gt, the caller's own claim, so a reference allele has no evaluation data and
 * every per-haplotype field reports "." for it. The evaluation lanes are indexed by calc_gt's
 * haplotypes, which matched_gt_is_swapped() reports may be the reverse of orig_gt's.
 * @param[in] out_fp Open file pointer to output VCF
 * @param[in] vi Variant index in this container
 * @param[in] sc_idx Supercluster index for SC field
 * @param[in] phase_block Phase block index for PB field
 * @param[in] phase_switch True if phase switched at this position
 * @param[in] phase_flip True if phase flipped (error) at this position
 * @param[in] query If true, format as query sample; if false, as truth sample
 */
void ctgVariants::print_var_sample(FILE* out_fp, int vi, int sc_idx, int phase_block,
        bool phase_switch, bool phase_flip, bool query /* = false */) {

    // a haploid record carries one bare allele; an unknown ploidy (0) is rendered as diploid
    int alleles = this->ploidies[vi] == 1 ? 1 : HAPS;
    const std::string gt = alleles == 1 ? "1" : gt_strs[this->orig_gts[vi]];

    bool swap = this->matched_gt_is_swapped(vi);
    std::string errtypes, credits, ref_eds, query_eds, match_types, sync_groups;
    for (int ai = 0; ai < alleles; ai++) {
        const std::string sep = ai ? "," : "";
        hap_t allele = hap_t(ai);

        // a reference allele was never evaluated, so it has no per-haplotype data to report
        if (!this->var_on_hap(vi, allele)) {
            errtypes += sep + "."; credits += sep + "."; ref_eds += sep + ".";
            query_eds += sep + "."; match_types += sep + "."; sync_groups += sep + ".";
            continue;
        }

        // get categorization
        hap_t hi = swap ? other_hap(allele) : allele;
        if (this->credit[hi][vi] == 1) {
            errtypes += sep + "TP"; match_types += sep + "gm";
        } else if (this->credit[hi][vi] == 0) {
            errtypes += sep + (query ? "FP" : "FN"); match_types += sep + ".";
        } else if (this->credit[hi][vi] >= g.credit_threshold) {
            errtypes += sep + "TP"; match_types += sep + "lm";
        } else {
            errtypes += sep + (query ? "FP" : "FN"); match_types += sep + "lm";
        }

        credits += sep + credit_str(this->credit[hi][vi]);
        ref_eds += sep + (this->ref_ed[hi][vi] == 0 ? "." :
                std::to_string(this->ref_ed[hi][vi]));
        query_eds += sep + (this->ref_ed[hi][vi] == 0 ? "." :
                std::to_string(this->query_ed[hi][vi]));
        sync_groups += sep + std::to_string(int(this->sync_group[hi][vi]));
    }

    fprintf(out_fp, "\t%s:%s:%s:%s:%s:%s:%d:%d:%s:%d:%d:%s:%s:%s:%s%s", gt.data(), errtypes.data(),
            credits.data(), ref_eds.data(), query_eds.data(), match_types.data(),
            int(this->var_quals[vi]), sc_idx, sync_groups.data(),
            this->phase_sets[vi], phase_block,
            query ? (phase_switch ? "1" : "0") : "." ,
            phase_strs[this->phases[vi]].data(),
            query ? (phase_flip ? "1" : "0") : "." ,
            ac_strs[this->ac_errtype[vi]].data(),
            query ? "\n" : "");
}

/**************************************************************************************************/

/**
 * @brief Constructs an empty variant data container defaulting to QUERY callset.
 */
variantData::variantData() : callset(QUERY) { ; }

/**
 * @brief Classifies a record's raw GT array into its parse-time genotype shape.
 *
 * Allele-index-agnostic: A and B stand for any alternate, so 1|2, 2|1, and 1|3 all classify as
 * compound heterozygous. Missing alleles are handled first, which leaves the diploid space
 * exhaustive over four cases on the two allele indices (both zero, exactly one zero, equal and
 * nonzero, distinct and nonzero) and the haploid space exhaustive over missing, zero, and nonzero.
 *
 * @param[in] gt Raw GT array as returned by bcf_get_format_int32(), unread when ngt is -1
 * @param[in] ngt Number of alleles in gt, or -1 when the record's VCF declares no GT tag
 * @return Parse-time genotype shape of the record
 * @throws ERROR Ploidy above 2, which the caller is expected to reject first
 */
gtparse_t classify_gt(const int32_t * gt, int ngt) {

    if (ngt == -1) return GT_PARSE_HAP_ALT; // no GT tag, assumed monoploid alternate
    if (ngt == 1) { // monoploid/haploid
        if (bcf_gt_is_missing(gt[0])) return GT_PARSE_HAP_MISSING;
        return bcf_gt_allele(gt[0]) ? GT_PARSE_HAP_ALT : GT_PARSE_HAP_REF;
    }
    if (ngt != 2) ERROR("classify_gt() expects monoploid/diploid GT, got ploidy %d", ngt);

    // distinguish a no-call (both alleles missing) from a half call (exactly one missing)
    const bool hap1_missing = bcf_gt_is_missing(gt[idx(HAP1)]);
    const bool hap2_missing = bcf_gt_is_missing(gt[idx(HAP2)]);
    if (hap1_missing && hap2_missing) return GT_PARSE_DIP_MISSING;
    if (hap1_missing || hap2_missing) return GT_PARSE_DIP_HALF_MISSING;

    const int allele1 = bcf_gt_allele(gt[idx(HAP1)]);
    const int allele2 = bcf_gt_allele(gt[idx(HAP2)]);
    if (!allele1 && !allele2) return GT_PARSE_DIP_HOM_REF;
    if (!allele1 || !allele2) return GT_PARSE_DIP_HET_ALT;
    if (allele1 == allele2)   return GT_PARSE_DIP_HOM_ALT;
    return GT_PARSE_DIP_CPD_HET_ALT;
}

/**
 * @brief Parses variants from a VCF file into a variantData container, with filtering and validation.
 * @param[in] vcf_fn Input VCF filename
 * @param[out] variant_data Container to populate with parsed variants
 * @param[in] reference Reference FASTA data for coordinate validation
 * @param[in] callset QUERY or TRUTH callset identifier
 * @throws ERROR A record htslib cannot parse, a record on a contig the header does not declare,
 *         a header declaring 'GT' with a type other than String, a header declaring other than one
 *         sample or a contig line without 'IDX' and 'length', an unsorted VCF, a variant of ploidy
 *         above 2, or a variant outside the reference contig
 * @throws WARNING Per-reason summary totals for records dropped or altered at parse time: no-call
 *         and half-call genotypes, records whose FORMAT column omits GT, unphased heterozygous
 *         genotypes, spanning deletions, reference calls, missing PS tags, records htslib parsed
 *         despite a non-critical error, oversized variants, overlapping variants, and complex
 *         variants split into INS + DEL
 */
void parse_variants(const std::string & vcf_fn,
        std::shared_ptr<variantData> variant_data,
        std::shared_ptr<fastaData> reference,
        callset_t callset) {

    // set reference fasta pointer
    variant_data->ref = reference;
    variant_data->filename = vcf_fn;

    variant_data->callset = callset;

    if (g.verbosity >= 1) INFO(" ");
    if (g.verbosity >= 1) INFO("%s[%s %d/%d] Parsing %s VCF%s '%s'", COLOR_PURPLE,
            callset == QUERY ? "Q" : "T", int(idx(TIME_READ)), int(idx(TIME_TOTAL))-1,
            callset_strs[callset].data(), 
            COLOR_WHITE, vcf_fn.data());
    htsFile* vcf = bcf_open(vcf_fn.data(), "r");

    // counters
    int nctg   = 0;                     // number of ctgs
    EnumArray<hap_t, EnumArray<edittype_t, int, EDITTYPE_SLOTS>, HAP_SLOTS> ntypes{};
    int n      = 0;                     // total number of records in file
    EnumArray<hap_t, int, HAP_SLOTS> npass{}; // records PASSing all filters

    // data
    bool print = g.verbosity >= 1;
    EnumArray<hap_t, int, HAP_SLOTS> prev_end =
            {{-g.cluster_min_gap*2, -g.cluster_min_gap*2}};
    EnumArray<hap_t, edittype_t, HAP_SLOTS> prev_type = {{TYPE_SUB, TYPE_SUB}};
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
    EnumArray<gtparse_t, int, GTPARSE_SLOTS> GT_counts{};
    int * gt      = NULL;
    bool gt_warn  = false;

    // phase set data for each call
    int PS_memsize   = 0;
    int nPS       = 0;
    int * PS      = NULL;
    bool PS_warn  = false;

    /* int gq_missing_total = 0; */
    int bcf_errcode_total = 0; // records htslib parsed despite a non-critical error
    int PS_missing_total = 0;
    int overlapping_var_total = 0;
    int spanning_del_total = 0;
    int no_gt_total = 0;          // records whose FORMAT omits a header-declared GT, dropped
    int unknown_allele_total = 0; // no-call records (.|. or .), dropped entirely
    int half_call_total = 0;      // half-call records (1|. or .|1), known allele kept
    int unphased_gt_total = 0;
    int too_large_var_total = 0;
    int multi_total = 0;
    int ref_call_total = 0;
    int complex_total = 0;
    int failed_filter_total = 0;
    
    // read header
    int read_ret = 0;      // bcf_read() return: 0 success, -1 end of file, < -1 critical error
    char errbuf[256] = ""; // bcf_strerror() decoding of rec->errcode
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
    
    while ((read_ret = bcf_read(vcf, hdr, rec)) != -1) { // -1 alone is end of file

        // a record htslib refuses to parse must not read as end of file, or a malformed callset is
        // scored as an empty one: a plausible recall of zero instead of a failure
        if (read_ret < -1) {
            bcf_strerror(rec->errcode, errbuf, sizeof(errbuf));
            ERROR("Failed to parse record %d of %s VCF '%s' at %s:%lld: %s", n+1,
                    callset_strs[callset].data(), vcf_fn.data(), bcf_seqname_safe(hdr, rec),
                    (long long)rec->pos, errbuf);
        }

        // htslib appends an undeclared contig to the header, so rid runs past ctgnames, which was
        // captured before the loop; indexing it would read out of bounds
        if (rec->rid < 0 || rec->rid >= nctg)
            ERROR("Contig '%s' in record %d of %s VCF '%s' is not declared in its header",
                    bcf_seqname_safe(hdr, rec), n+1, callset_strs[callset].data(),
                    vcf_fn.data());

        // every other errcode is recoverable, and htslib parsed the rest of the record; aborting
        // would reject the many real VCFs that carry an undeclared INFO or FORMAT tag
        if (rec->errcode) {
            if (g.verbosity > 1) {
                bcf_strerror(rec->errcode, errbuf, sizeof(errbuf));
                WARN("Record %d of %s VCF parsed with an error at %s:%lld: %s", n+1,
                        callset_strs[callset].data(), bcf_seqname_safe(hdr, rec),
                        (long long)rec->pos, errbuf);
            }
            bcf_errcode_total++;
        }

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
                prev_end = {{-g.cluster_min_gap*2, -g.cluster_min_gap*2}};
                prev_type = {{TYPE_SUB, TYPE_SUB}};
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
        // the spec fixes GT's type as String, so no retry is possible the way GQ retries as float
        } else if (ngt == -2) { // header declares GT with a type the spec forbids
            ERROR("%s VCF header declares 'GT' with a type other than String, at %s:%lld",
                    callset_strs[callset].data(), ctg.data(), (long long)rec->pos);

        // legal VCF: FORMAT is per-record and GT is not mandatory. No genotype can be assumed
        // without fabricating an allele that would then be scored as a TP or FN
        } else if (ngt == -3) { // GT declared in the header but absent from this record's FORMAT
            if (g.verbosity > 1)
                WARN("No GT tag in %s VCF at %s:%lld, skipping",
                        callset_strs[callset].data(), ctg.data(), (long long)rec->pos);
            no_gt_total++;
            continue;

        } else if (ngt <= 0) { // -4 and 0, neither reachable: bcf_read() rejects the record first
            ERROR("Failed to read %s GT at %s:%lld (htslib returned %d)",
                    callset_strs[callset].data(), ctg.data(), (long long)rec->pos, ngt);
        }

        // record this record's ploidy; mixed ploidy within a contig is legitimate, as on a chrX
        // carrying both PAR (diploid) and non-PAR (haploid) calls, so nothing is enforced here
        int ctg_idx = std::find(variant_data->contigs.begin(), variant_data->contigs.end(), ctg)
                - variant_data->contigs.begin();
        variant_data->observed_ploidies[ctg_idx].insert(std::abs(ngt));

        // parse genotype info
        if (ngt > 2) // polyploid, rejected before classify_gt() sees it
            ERROR("Expected monoploid/diploid %s VCF, found variant with ploidy %d",
                    callset_strs[callset].data(), ngt);
        const gtparse_t parse_gt = classify_gt(gt, ngt);
        GT_counts[parse_gt]++;

        // both alleles equal allows setting N/N to 1/1 later
        const bool same = parse_gt == GT_PARSE_DIP_HOM_REF || parse_gt == GT_PARSE_DIP_HOM_ALT;

        // count missing alleles once per record, not once per haplotype
        if (parse_gt == GT_PARSE_DIP_MISSING || parse_gt == GT_PARSE_HAP_MISSING) {
            if (g.verbosity > 1)
                WARN("Variant with no known alleles (.|.) in %s VCF at %s:%lld, skipping",
                    callset_strs[callset].data(), ctg.data(), (long long)rec->pos);
            unknown_allele_total += 1;
        } else if (parse_gt == GT_PARSE_DIP_HALF_MISSING) {
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
        EnumArray<hap_t, int, HAP_SLOTS> rec_prev_end = prev_end;
        EnumArray<hap_t, edittype_t, HAP_SLOTS> rec_prev_type = prev_type;

        // parse variant type
        for (int hi = 0; hi < std::abs(ngt); hi++) { // allow single-allele chrX, chrY
            hap_t hap = static_cast<hap_t>(hi);

            // set simplified GT (0|1, 1|0, or 1|1), (0|0 and .|. skipped later)
            gt_t simple_gt = hap == HAP2 ? GT_REF_ALT : GT_ALT_REF; // 0|1 or 1|0 default
            if (same) simple_gt = GT_ALT_ALT; // overwrite 1|1 if both agree

            // get ref and allele, skipping ref query
            std::string ref = rec->d.allele[0];
            int alt_idx = ngt < 0 ? 1 : bcf_gt_allele(gt[idx(hap)]); // if no GT, assume 1
            if (alt_idx < 0) continue; // missing allele (.), counted once per record above
            if (alt_idx == 0) continue; // nothing to do if reference
            std::string alt = rec->d.allele[alt_idx];

            // uppercase before any comparison, since soft-masked reference sequence reaches us as
            // lowercase and case must not decide whether alleles match
            std::transform(ref.begin(), ref.end(), ref.begin(), ::toupper);
            std::transform(alt.begin(), alt.end(), alt.begin(), ::toupper);

            // skip unphased heterozygous variants (1/1 is allowed, 0/1 is not)
            if (ngt == 2 && !same && !bcf_gt_is_phased(gt[idx(HAP2)])) { // only HAP2 is set, not sure why...
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
            edittype_t type = TYPE_REF;
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
            if (simple_gt == GT_ALT_ALT && (rec_prev_end[other_hap(hap)] > pos ||
                    (rec_prev_end[other_hap(hap)] == pos && rec_prev_type[other_hap(hap)] == TYPE_INS &&
                     type == TYPE_INS))) {
                simple_gt = hap == HAP2 ? GT_REF_ALT : GT_ALT_REF;
            }

            // add to haplotype-specific query info
            int rec_idx = n - 1; // 0-based ordinal of this record within the input VCF
            uint8_t ploidy = uint8_t(std::abs(ngt));
            // both CPX halves derive from the same original allele, so they share alt_idx
            if (type == TYPE_CPX) { // split CPX into INS+DEL
                variant_data->variants[hap][ctg]->add_var(var_fields{.pos = pos, .rlen = 0, // INS
                    .type = TYPE_INS, .loc = loc, .ref = "", .alt = alt,
                    .orig_gt = simple_gt, .gt_qual = float(ngq ? gq[0]:0),
                    .var_qual = vq, .phase_set = phase_set,
                    .rec_idx = rec_idx, .alt_idx = alt_idx, .ploidy = ploidy});
                variant_data->variants[hap][ctg]->add_var(var_fields{.pos = pos, .rlen = rlen, // DEL
                    .type = TYPE_DEL, .loc = loc, .ref = ref, .alt = "",
                    .orig_gt = simple_gt, .gt_qual = float(ngq ? gq[0]:0),
                    .var_qual = vq, .phase_set = phase_set,
                    .rec_idx = rec_idx, .alt_idx = alt_idx, .ploidy = ploidy});
                complex_total++;
            } else {
                variant_data->variants[hap][ctg]->add_var(var_fields{.pos = pos, .rlen = rlen,
                        .type = type, .loc = loc, .ref = ref, .alt = alt,
                        .orig_gt = simple_gt, .gt_qual = float(ngq ? gq[0]:0),
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
    for (gtparse_t gt : EnumRange<gtparse_t, GTPARSE_SLOTS>{}) {
        if (print && GT_counts[gt]) INFO("    %3s: %i", gtparse_strs[gt].data(), GT_counts[gt]);
    }
    if (print) INFO(" ");

    if (bcf_errcode_total)
        WARN("%d records with parse errors in %s VCF, kept",
            bcf_errcode_total, callset_strs[callset].data());

    if (PS_missing_total) 
        WARN("%d variants missing PS tags in %s VCF, kept",
            PS_missing_total, callset_strs[callset].data());

    multi_total = GT_counts[GT_PARSE_DIP_HOM_ALT] + GT_counts[GT_PARSE_DIP_CPD_HET_ALT];
    if (multi_total && print)
        INFO("%d homozygous and multi-allelic variants in %s VCF, split for evaluation",
            multi_total, callset_strs[callset].data());

    if (no_gt_total)
        WARN("%d variants with no GT field in %s VCF, skipped",
            no_gt_total, callset_strs[callset].data());

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
        for (hap_t h : EnumRange<hap_t, HAP_SLOTS>{}) {
            if (print) INFO("    Haplotype %i", int(idx(h))+1);
            for (edittype_t t : EnumRange<edittype_t, EDITTYPE_SLOTS>{}) {
                if (print) INFO("      %s: %i", type_strs[t].data(), ntypes[h][t]);
            }
        }
        if (print) INFO(" ");
    } else { // summarize
        for (edittype_t t : EnumRange<edittype_t, EDITTYPE_SLOTS>{}) {
            if (print && ntypes[HAP1][t] + ntypes[HAP2][t])
                INFO("    %s: %i", type_strs[t].data(), ntypes[HAP1][t] + ntypes[HAP2][t]);
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
