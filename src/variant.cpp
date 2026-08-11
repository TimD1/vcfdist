/**
 * @file variant.cpp
 * @brief Per-contig and per-callset variant containers with VCF parsing and output utilities.
 */
#include <algorithm>
#include <chrono>
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


/* Match tier criteria ****************************************************************************/


/**
 * @brief Maps an allele count error type onto the query's allele count relative to the truth's.
 *
 * A total mapping, so that no allele count error type falls outside the match tier ladder. The
 * 0 -> N rows are pure query false positives and the N -> 0 rows pure truth false negatives; both
 * are counted as errors in a direction rather than left undefined, since only ALLELE_COUNT_EQUAL
 * reaches the gm tier and the direction is what separates an allele error from a missed call.
 *
 * A 0 -> N row does not imply zero credit. A query heterozygote whose best haplotype scored below
 * --credit-threshold has no matched haplotype, so it lands on AC_ERR_0_TO_1 while its
 * MaxAlleleCredit is still CREDIT_NONZERO; it reaches the lm tier and stops at am.
 *
 * @param[in] ac_errtype Allele count error type, from ctgVariants::ac_errtype
 * @return Query allele count relative to the truth allele count
 * @throws ERROR if the allele count error type is AC_UNKNOWN
 */
allelecount_t ac_errtype_to_allele_count(ac_errtype_t ac_errtype) {
    switch (ac_errtype) {
        case AC_ERR_0_TO_1:
        case AC_ERR_0_TO_2:
        case AC_ERR_1_TO_2: return ALLELE_COUNT_GAIN;
        case AC_ERR_1_TO_1:
        case AC_ERR_2_TO_2: return ALLELE_COUNT_EQUAL;
        case AC_ERR_1_TO_0:
        case AC_ERR_2_TO_0:
        case AC_ERR_2_TO_1: return ALLELE_COUNT_LOSS;
        case AC_UNKNOWN: break;
    }
    // fix_allele_counts() already errors on AC_UNKNOWN for both callsets, so this is unreachable
    ERROR("Unknown allele count error type in ac_errtype_to_allele_count()");
}


/**
 * @brief Returns the most stringent match tier the three criteria jointly satisfy.
 *
 * Each tier adds one criterion to the one below it, so the result is monotone by construction:
 * every tier's conditions are a superset of the conditions of every tier beneath it.
 *
 * MaxAlleleCredit supplies both lower rungs and MinAlleleCredit neither. Since Min <= Max always,
 * Min is the stronger predicate, so putting it on the looser rung would invert the ladder: a query
 * 1|1 against a truth 0|1 (Max = PASS, Min = ZERO) would satisfy am while failing lm, ranking a
 * genotype error below a weak partial match. Two thresholds on Max are monotone instead, and
 * discriminate on heterozygotes too, where the single carried haplotype makes Min == Max.
 *
 * @param[in] max_credit Highest per-haplotype credit, bucketed against --credit-threshold
 * @param[in] allele_count Query allele count relative to the truth allele count
 * @param[in] phase_match Whether the variant's phasing matches the phasing its block chose
 * @return Most stringent tier satisfied, or MATCH_NONE if not even lm is reached
 */
matchtier_t match_tier(credit_t max_credit, allelecount_t allele_count, phasematch_t phase_match) {
    if (max_credit == CREDIT_ZERO) return MATCH_NONE;
    if (max_credit != CREDIT_PASS) return MATCH_LM;
    if (allele_count != ALLELE_COUNT_EQUAL) return MATCH_AM;
    if (phase_match != PHASEMATCH_CORRECT && phase_match != PHASEMATCH_NOT_HETEROZYGOUS)
        return MATCH_GM;
    return MATCH_PM;
}


/**
 * @brief Returns a variant's highest per-haplotype credit, bucketed against --credit-threshold.
 *
 * The maximum runs over both haplotype lanes rather than only the carried ones, which needs no
 * special case for a heterozygote: the lane it does not carry was never evaluated and holds zero.
 *
 * @param[in] vi Variant index
 * @return CREDIT_ZERO, CREDIT_NONZERO, or CREDIT_PASS
 */
credit_t ctgVariants::get_max_allele_credit(int vi) const {
    float max_credit = 0;
    for (hap_t hap : EnumRange<hap_t, HAP_SLOTS>{}) {
        max_credit = std::max(max_credit, this->credit[hap][vi]);
    }
    if (max_credit <= 0) return CREDIT_ZERO;
    // --credit-threshold is validated into (0, 1], so zero credit can never reach CREDIT_PASS
    return max_credit >= g.credit_threshold ? CREDIT_PASS : CREDIT_NONZERO;
}


/**
 * @brief Returns whether a variant's alignment phasing matches the phasing its block chose.
 *
 * The homozygous reference, the homozygous alternate, and the haploid call are all
 * PHASEMATCH_NOT_HETEROZYGOUS, since none of them occupies a distinguishable pair of haplotypes.
 * The homozygous reference cannot currently reach here, as parse_variants() stores one variant per
 * non-reference allele, but it is named rather than left to fall through. Naming the homozygous
 * genotypes rather than the heterozygous ones is deliberate: a heterozygous genotype added to gt_t
 * later falls through to the phase comparison, which is the answer it wants.
 *
 * The criterion is query-side: phases and pb_phases are populated during phasing for query
 * variants only, so a truth heterozygote reports PHASEMATCH_UNPHASED. That costs nothing, because
 * a truth-only site has no query record whose phase could be verified and fails at the lm rung
 * regardless.
 *
 * @param[in] vi Variant index
 * @return PHASEMATCH_CORRECT, PHASEMATCH_INCORRECT, PHASEMATCH_UNPHASED, or
 *         PHASEMATCH_NOT_HETEROZYGOUS
 */
phasematch_t ctgVariants::get_phase_match(int vi) const {
    // homozygous and haploid variants keep the PHASE_NONE default that add_var() sets, so they must
    // be recognized before it, or every one of them would be reported as merely unphased
    if (this->orig_gts[vi] == GT_REF_REF || this->orig_gts[vi] == GT_ALT_ALT ||
            this->ploidies[vi] == PLOIDY_HAPLOID) {
        return PHASEMATCH_NOT_HETEROZYGOUS;
    }
    if (this->phases[vi] == PHASE_NONE || this->pb_phases[vi] == PHASE_NONE) {
        return PHASEMATCH_UNPHASED;
    }
    return this->phases[vi] == this->pb_phases[vi] ? PHASEMATCH_CORRECT : PHASEMATCH_INCORRECT;
}


/**
 * @brief Returns the most stringent match tier a variant satisfies.
 * @param[in] vi Variant index
 * @return Most stringent tier satisfied, or MATCH_NONE if not even lm is reached
 * @throws ERROR if the variant's allele count error type is AC_UNKNOWN
 */
matchtier_t ctgVariants::get_match_tier(int vi) const {
    return match_tier(this->get_max_allele_credit(vi),
            ac_errtype_to_allele_count(this->ac_errtype[vi]),
            this->get_phase_match(vi));
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
 * @brief Builds the summary VCF header, declaring every FORMAT field and the TRUTH/QUERY samples.
 * @param[in] contigs Contig names, in the order records are written
 * @param[in] lengths Contig lengths, parallel to contigs
 * @return Header owning its own memory, to be released by the caller with bcf_hdr_destroy()
 * @throws ERROR The header cannot be allocated, a header line htslib rejects, a sample htslib
 *         rejects, or a header htslib cannot synchronize
 */
bcf_hdr_t* summary_vcf_header(const std::vector<std::string> & contigs,
        const std::vector<int> & lengths) {

    // bcf_hdr_init() supplies the ##fileformat line
    bcf_hdr_t* hdr = bcf_hdr_init("w");
    if (hdr == NULL) ERROR("Failed to allocate summary VCF header");

    const std::chrono::time_point<std::chrono::system_clock> now{std::chrono::system_clock::now()};
    time_t tt = std::chrono::system_clock::to_time_t(now);
    tm local_time = *localtime(&tt);
    char file_date[32];
    snprintf(file_date, sizeof(file_date), "##fileDate=%04d%02d%02d", local_time.tm_year + 1900,
            local_time.tm_mon + 1, local_time.tm_mday);

    // The per-haplotype fields carry one value per allele of the sample's GT, which is what VCF
    // 4.4's Number=P declares. BCF_VL_P only reaches htslib in 1.23, so a consumer on any older
    // bcftools or pysam would report a cardinality error; Number=. produces identical records and
    // merely gives up the declared cardinality, so the count and order are stated here instead.
    const std::string per_allele = " One value per allele of this sample's GT, in GT allele order, "
            "'.' for a reference allele.";
    std::vector<std::string> lines = {file_date, "##CL=" + g.cmd};
    for (size_t i = 0; i < contigs.size(); i++) {
        lines.push_back("##contig=<ID=" + contigs[i] + ",length=" +
                std::to_string(lengths[i]) + ">");
    }
    // every record PASSes, and htslib rejects a filter its header does not declare; bcf_hdr_init()
    // declares PASS itself, and appending an ID it already holds is a no-op rather than a duplicate
    lines.push_back("##FILTER=<ID=PASS,Description=\"All filters passed\">");
    lines.push_back("##FORMAT=<ID=GT,Number=1,Type=String,Description=\"GenoType\">");
    lines.push_back("##FORMAT=<ID=BD,Number=.,Type=String,Description=\"Benchmark Decision for call (TP/FP/FN)." + per_allele + "\">");
    lines.push_back("##FORMAT=<ID=BC,Number=.,Type=Float,Description=\"Benchmark Credit (on the interval [0,1], based on sync group edit distance)." + per_allele + "\">");
    lines.push_back("##FORMAT=<ID=RD,Number=.,Type=Integer,Description=\"Reference edit Distance from truth within current sync group." + per_allele + "\">");
    lines.push_back("##FORMAT=<ID=QD,Number=.,Type=Integer,Description=\"Query edit Distance from truth within current sync group." + per_allele + "\">");
    lines.push_back("##FORMAT=<ID=BK,Number=.,Type=String,Description=\"BenchmarK category ('gm' if credit == 1, 'lm' if credit > 0, else '.')." + per_allele + "\">");
    lines.push_back("##FORMAT=<ID=QQ,Number=1,Type=Float,Description=\"variant Quality\">");
    lines.push_back("##FORMAT=<ID=SC,Number=1,Type=Integer,Description=\"SuperCluster (index in contig)\">");
    lines.push_back("##FORMAT=<ID=SG,Number=.,Type=Integer,Description=\"Sync Group (index in supercluster, for credit assignment)." + per_allele + "\">");
    lines.push_back("##FORMAT=<ID=PS,Number=1,Type=Integer,Description=\"Phase Set identifier (input, per-variant)\">");
    lines.push_back("##FORMAT=<ID=PB,Number=1,Type=Integer,Description=\"Phase Block (output, per-supercluster, index in contig)\">");
    lines.push_back("##FORMAT=<ID=BS,Number=1,Type=Integer,Description=\"Block Phase: 0 = PHASE_KEEP, 1 = PHASE_SWAP)\">");
    lines.push_back("##FORMAT=<ID=VP,Number=1,Type=Integer,Description=\"Variant Phase: 0 = PHASE_ORIG, 1 = PHASE_SWAP, . = PHASE_NONE)\">");
    lines.push_back("##FORMAT=<ID=FE,Number=1,Type=Integer,Description=\"Flip Error (a per-supercluster error)\">");
    lines.push_back("##FORMAT=<ID=GE,Number=1,Type=String,Description=\"Genotype Error ('+' if 0/1 truth -> 1/1 query, '-' if 1/1 truth -> 0/1 query, '.' otherwise)\">");
    for (const std::string & line : lines) {
        if (bcf_hdr_append(hdr, line.data()) != 0)
            ERROR("Failed to add summary VCF header line '%s'", line.data());
    }

    // the samples are added in the order their columns are written
    for (const char* sample : {"TRUTH", "QUERY"}) {
        if (bcf_hdr_add_sample(hdr, sample) < 0)
            ERROR("Failed to add sample '%s' to summary VCF header", sample);
    }
    if (bcf_hdr_sync(hdr) < 0) ERROR("Failed to synchronize summary VCF header");
    return hdr;
}


/**
 * @brief Sets the fixed VCF fields (CHROM, POS, ID, REF, ALT, QUAL, FILTER) of one record.
 * @param[in] hdr Summary VCF header, which must declare this contig
 * @param[in,out] rec Cleared record to fill
 * @param[in] ref Reference FASTA data for retrieving flanking bases for indels
 * @param[in] ctg Contig name
 * @param[in] idx Variant index in this container
 * @throws ERROR The contig is not declared in the header
 * @throws ERROR An INS/DEL sits at the contig start (0-based pos 0), leaving no preceding base to anchor
 * @throws ERROR The variant type is not TYPE_SUB, TYPE_INS, or TYPE_DEL
 * @throws ERROR htslib rejects the record's FILTER or alleles
 */
void ctgVariants::set_var_record(const bcf_hdr_t* hdr, bcf1_t* rec,
        std::shared_ptr<fastaData> ref, const std::string & ctg, int idx) const {

    rec->rid = bcf_hdr_name2id(hdr, ctg.data());
    if (rec->rid < 0) ERROR("Contig '%s' is not declared in the summary VCF header", ctg.data());
    bcf_float_set_missing(rec->qual);
    if (bcf_add_filter(hdr, rec, bcf_hdr_id2int(hdr, BCF_DT_ID, "PASS")) < 0)
        ERROR("Failed to set FILTER on summary VCF record at %s:%d", ctg.data(), this->poss[idx]);

    std::string ref_allele, alt_allele;
    switch (this->types[idx]) {
    case TYPE_SUB:
        rec->pos = this->poss[idx];
        ref_allele = this->refs[idx];
        alt_allele = this->alts[idx];
        break;
    case TYPE_INS:
    case TYPE_DEL: {
        // INS/DEL are left-anchored on the preceding reference base; at contig start (0-based
        // pos 0) there is no preceding base, so guard against the out-of-bounds read of index -1
        if (this->poss[idx] == 0)
            ERROR("Cannot left-anchor INS/DEL at contig start (0-based pos 0) on '%s' in set_var_record",
                    ctg.data());
        rec->pos = this->poss[idx] - 1;
        char ref_base = ref->fasta.at(ctg)[this->poss[idx]-1];
        ref_allele = ref_base + this->refs[idx];
        alt_allele = ref_base + this->alts[idx];
        break;
    }
    default:
        ERROR("set_var_record not implemented for type %d", static_cast<int>(this->types[idx]));
    }

    const char* alleles[2] = {ref_allele.data(), alt_allele.data()};
    if (bcf_update_alleles(hdr, rec, alleles, 2) < 0)
        ERROR("Failed to set alleles on summary VCF record at %s:%d", ctg.data(), this->poss[idx]);
}


/**
 * @brief Encodes the GT a sample reports for one variant, as the caller itself genotyped it.
 * @param[in] orig_gt The caller's own genotype, never vcfdist's recovered matched_gt
 * @param[in] ploidy Variant ploidy
 * @return Phased allele indices: one for a haploid call, otherwise the diploid pair
 */
static std::vector<int32_t> genotype_alleles(gt_t orig_gt, ploidy_t ploidy) {
    if (ploidy == PLOIDY_HAPLOID) return {bcf_gt_phased(1)};
    return {bcf_gt_phased(orig_gt == GT_ALT_REF || orig_gt == GT_ALT_ALT ? 1 : 0),
            bcf_gt_phased(orig_gt == GT_REF_ALT || orig_gt == GT_ALT_ALT ? 1 : 0)};
}


/**
 * @brief Returns the FORMAT values of a sample that made no call at a locus.
 * @param[in] sc_idx Supercluster index for the SC field
 * @param[in] phase_block Phase block index for the PB field
 * @return Values reporting '.' for every field but SC and PB, which are locus-wide
 */
sample_fields empty_sample_fields(int sc_idx, int phase_block) {
    sample_fields fields;
    bcf_float_set_missing(fields.qq);
    fields.sc = sc_idx;
    fields.ps = bcf_int32_missing;
    fields.pb = phase_block;
    fields.bs = bcf_int32_missing;
    fields.vp = bcf_int32_missing;
    fields.fe = bcf_int32_missing;
    return fields;
}


/**
 * @brief Returns one sample's FORMAT values for a variant it called.
 *
 * One record is written per variant rather than per haplotype, so the per-haplotype fields (BD, BC,
 * RD, QD, BK, SG) hold one value per haplotype: two for a diploid record, one for a haploid one. GT
 * reports orig_gt, the caller's own claim, so a haplotype carrying the reference allele has no
 * evaluation data and every per-haplotype field reports '.' for it. The evaluation lanes are keyed
 * by matched_gt's haplotypes, which matched_gt_is_swapped() reports may be the reverse of orig_gt's.
 * @param[in] vi Variant index in this container
 * @param[in] sc_idx Supercluster index for the SC field
 * @param[in] phase_block Phase block index for the PB field
 * @param[in] phase_switch True if phase switched at this position
 * @param[in] phase_flip True if phase flipped (error) at this position
 * @param[in] query If true, report as the query sample; if false, as the truth sample
 * @return This sample's FORMAT values, htslib-encoded
 */
sample_fields ctgVariants::var_sample_fields(int vi, int sc_idx, int phase_block,
        bool phase_switch, bool phase_flip, bool query /* = false */) const {

    // ploidy is the count of genotype alleles, so it is also how many haplotypes to report on
    ploidy_t ploidy = this->ploidies[vi];
    int haps = int(idx(ploidy));

    sample_fields fields;
    fields.gt = genotype_alleles(this->orig_gts[vi], ploidy);
    // QQ was printed with %d before this file wrote records through htslib, so it stays truncated
    fields.qq = float(int(this->var_quals[vi]));
    fields.sc = sc_idx;
    fields.ps = this->phase_sets[vi];
    fields.pb = phase_block;
    fields.bs = query ? (phase_switch ? 1 : 0) : bcf_int32_missing;
    fields.vp = this->phases[vi] == PHASE_NONE ?
            bcf_int32_missing : int32_t(idx(this->phases[vi]));
    fields.fe = query ? (phase_flip ? 1 : 0) : bcf_int32_missing;
    fields.ge = ac_strs[this->ac_errtype[vi]];

    bool swap = this->matched_gt_is_swapped(vi);
    float missing_credit;
    bcf_float_set_missing(missing_credit);
    std::string errtypes, match_types;
    for (int hap_idx = 0; hap_idx < haps; hap_idx++) {
        const std::string sep = hap_idx ? "," : "";
        hap_t hi = hap_t(hap_idx);

        // this haplotype carries the reference allele, so it was never evaluated
        if (!this->var_on_hap(vi, hi)) {
            errtypes += sep + "."; match_types += sep + ".";
            fields.bc.push_back(missing_credit);
            fields.rd.push_back(bcf_int32_missing);
            fields.qd.push_back(bcf_int32_missing);
            fields.sg.push_back(bcf_int32_missing);
            continue;
        }

        // the evaluation lanes are keyed by matched_gt's haplotypes, not orig_gt's
        hap_t hi_matched = swap ? other_hap(hi) : hi;

        // get categorization
        if (this->credit[hi_matched][vi] == 1) {
            errtypes += sep + "TP"; match_types += sep + "gm";
        } else if (this->credit[hi_matched][vi] == 0) {
            errtypes += sep + (query ? "FP" : "FN"); match_types += sep + ".";
        } else if (this->credit[hi_matched][vi] >= g.credit_threshold) {
            errtypes += sep + "TP"; match_types += sep + "lm";
        } else {
            errtypes += sep + (query ? "FP" : "FN"); match_types += sep + "lm";
        }

        fields.bc.push_back(this->credit[hi_matched][vi]);
        fields.rd.push_back(this->ref_ed[hi_matched][vi] == 0 ?
                bcf_int32_missing : this->ref_ed[hi_matched][vi]);
        fields.qd.push_back(this->ref_ed[hi_matched][vi] == 0 ?
                bcf_int32_missing : this->query_ed[hi_matched][vi]);
        fields.sg.push_back(this->sync_group[hi_matched][vi]);
    }
    fields.bd = errtypes;
    fields.bk = match_types;
    return fields;
}


/** @brief Sets one value to the missing value of its type. */
static void set_missing(int32_t & value) { value = bcf_int32_missing; }
static void set_missing(float & value) { bcf_float_set_missing(value); }

/** @brief Sets one value to the end-of-vector marker of its type. */
static void set_vector_end(int32_t & value) { value = bcf_int32_vector_end; }
static void set_vector_end(float & value) { bcf_float_set_vector_end(value); }

/**
 * @brief Concatenates two samples' per-allele values into one buffer of a shared value count.
 *
 * htslib stores the same number of values for every sample, so a haploid sample beside a diploid
 * one is padded with the end-of-vector marker, which the writer prints as a shorter list. A sample
 * with no values at all still occupies one slot, holding the missing value.
 * @param[in] truth The truth sample's values
 * @param[in] query The query sample's values
 * @return The truth sample's padded values followed by the query sample's
 */
template <typename T>
static std::vector<T> pad_per_allele(const std::vector<T> & truth, const std::vector<T> & query) {
    const std::vector<T>* samples[2] = {&truth, &query};
    size_t n = std::max(std::max(truth.size(), query.size()), size_t(1));
    std::vector<T> values(2 * n);
    for (size_t si = 0; si < 2; si++) {
        for (size_t i = 0; i < n; i++) {
            T & value = values[si*n + i];
            if (i < samples[si]->size()) value = (*samples[si])[i];
            else if (i == 0) set_missing(value);
            else set_vector_end(value);
        }
    }
    return values;
}

/** @brief Sets one integer FORMAT field of a record, holding one value per sample per allele. */
static void update_format(const bcf_hdr_t* hdr, bcf1_t* rec, const char* key,
        const std::vector<int32_t> & values) {
    if (bcf_update_format_int32(hdr, rec, key, values.data(), int(values.size())) < 0)
        ERROR("Failed to set FORMAT/%s on summary VCF record", key);
}

/** @brief Sets one float FORMAT field of a record, holding one value per sample per allele. */
static void update_format(const bcf_hdr_t* hdr, bcf1_t* rec, const char* key,
        const std::vector<float> & values) {
    if (bcf_update_format_float(hdr, rec, key, values.data(), int(values.size())) < 0)
        ERROR("Failed to set FORMAT/%s on summary VCF record", key);
}

/** @brief Sets one string FORMAT field of a record, holding one string per sample. */
static void update_format(const bcf_hdr_t* hdr, bcf1_t* rec, const char* key,
        const std::string & truth, const std::string & query) {
    const char* values[2] = {truth.data(), query.data()};
    if (bcf_update_format_string(hdr, rec, key, values, 2) < 0)
        ERROR("Failed to set FORMAT/%s on summary VCF record", key);
}


/**
 * @brief Sets every FORMAT field of one record from the two samples' values.
 *
 * Fields are added in FORMAT declaration order, which is the order htslib writes them in.
 * @param[in] hdr Summary VCF header, which must declare every field
 * @param[in,out] rec Record whose fixed fields are already set
 * @param[in] truth The truth sample's values
 * @param[in] query The query sample's values
 * @throws ERROR htslib rejects any field's values
 */
void set_record_samples(const bcf_hdr_t* hdr, bcf1_t* rec,
        const sample_fields & truth, const sample_fields & query) {

    // a sample with no call reports one missing allele, the '.' genotype
    std::vector<int32_t> truth_gt = truth.gt, query_gt = query.gt;
    if (truth_gt.empty()) truth_gt.push_back(bcf_gt_missing);
    if (query_gt.empty()) query_gt.push_back(bcf_gt_missing);
    const std::vector<int32_t> gt = pad_per_allele(truth_gt, query_gt);
    if (bcf_update_genotypes(hdr, rec, gt.data(), int(gt.size())) < 0)
        ERROR("Failed to set FORMAT/GT on summary VCF record");

    update_format(hdr, rec, "BD", truth.bd, query.bd);
    update_format(hdr, rec, "BC", pad_per_allele(truth.bc, query.bc));
    update_format(hdr, rec, "RD", pad_per_allele(truth.rd, query.rd));
    update_format(hdr, rec, "QD", pad_per_allele(truth.qd, query.qd));
    update_format(hdr, rec, "BK", truth.bk, query.bk);
    update_format(hdr, rec, "QQ", std::vector<float>{truth.qq, query.qq});
    update_format(hdr, rec, "SC", std::vector<int32_t>{truth.sc, query.sc});
    update_format(hdr, rec, "SG", pad_per_allele(truth.sg, query.sg));
    update_format(hdr, rec, "PS", std::vector<int32_t>{truth.ps, query.ps});
    update_format(hdr, rec, "PB", std::vector<int32_t>{truth.pb, query.pb});
    update_format(hdr, rec, "BS", std::vector<int32_t>{truth.bs, query.bs});
    update_format(hdr, rec, "VP", std::vector<int32_t>{truth.vp, query.vp});
    update_format(hdr, rec, "FE", std::vector<int32_t>{truth.fe, query.fe});
    update_format(hdr, rec, "GE", truth.ge, query.ge);
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
    hts_pos_t prev_pos = -1; // previous record's position within the current contig
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
                prev_pos = -1;
            }
        }

        // a record moving backwards within a contig is otherwise dropped silently by the overlap
        // filter below, leaving a run that succeeds with every denominator quietly wrong
        if (rec->pos < prev_pos)
            ERROR("Unsorted %s VCF '%s', record %d at %s:%lld precedes position %lld",
                    callset_strs[callset].data(), vcf_fn.data(), n+1, ctg.data(),
                    (long long)rec->pos, (long long)prev_pos);
        prev_pos = rec->pos;

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

            // check that variant (original representation) is in region of interest; with no BED
            // supplied there are no regions to be outside of, so every variant is evaluated
            bedloc_t loc = g.bed_exists ?
                    g.bed.contains(ctg, rec->pos, rec->pos + reflen, type) : BED_INSIDE;
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
            // ngt is 1 or 2 by here: a polyploid record errored above, and a record whose VCF
            // declares no GT tag reports -1 and is treated as the monoploid call it is assumed to be
            ploidy_t ploidy = ploidy_t(std::abs(ngt));
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
