/**
 * @file variant.cpp
 * @brief Per-contig and per-callset variant containers with VCF parsing and output utilities.
 */
#include <algorithm>
#include <chrono>
#include <set>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <vector>
#include <cmath>
#include <cstdio>
#include <cstdlib>

#include "htslib/vcf.h"

#include "variant.h"
#include "print.h"
#include "dist.h"

/* Source column retention ************************************************************************/

/** @brief 0-based column indices of a record in a single-sample VCF. */
enum vcf_col { CHROM_COL, POS_COL, ID_COL, REF_COL, ALT_COL, QUAL_COL, FILTER_COL, INFO_COL,
        FORMAT_COL, SAMPLE_COL, VCF_COLS };

/**
 * @brief Splits a string on a delimiter, keeping empty fields.
 * @param[in] str String to split
 * @param[in] delim Delimiter to split on
 * @return Every field between delimiters, in order
 */
static std::vector<std::string> split(const std::string & str, char delim) {
    std::vector<std::string> fields;
    size_t start = 0;
    for (size_t i = 0; i <= str.size(); i++) {
        if (i == str.size() || str[i] == delim) {
            fields.push_back(str.substr(start, i - start));
            start = i + 1;
        }
    }
    return fields;
}

/** @brief FORMAT keys the summary VCF writer emits itself, in the order it emits them. */
static const std::string FIXED_FMT_KEYS = "GT:BD:BC:RD:QD:BK:QQ:SC:SG:PS:PB:BS:VP:FE:GE";

/**
 * @brief Splits the fixed FORMAT key list into a set, so source keys can be checked against it.
 * @return One entry per key the summary VCF writer emits itself
 */
static std::set<std::string> fixed_fmt_key_set() {
    std::vector<std::string> keys = split(FIXED_FMT_KEYS, ':');
    return std::set<std::string>(keys.begin(), keys.end());
}

/** @brief Derived from FIXED_FMT_KEYS so the two can never disagree about what the writer emits. */
static const std::set<std::string> FIXED_FMT_KEY_SET = fixed_fmt_key_set();

/**
 * @brief Returns a field's declared length class, or BCF_VL_FIXED when it is undeclared.
 * @param[in] hdr Header of the VCF the field was declared in
 * @param[in] line_type Header line type the field is declared on (BCF_HL_INFO or BCF_HL_FMT)
 * @param[in] key Field name
 * @return The htslib length class (BCF_VL_FIXED, BCF_VL_VAR, BCF_VL_A, BCF_VL_R, or BCF_VL_G)
 */
static int length_class(const bcf_hdr_t * hdr, int line_type, const std::string & key) {
    int id = bcf_hdr_id2int(hdr, BCF_DT_ID, key.data());
    if (id < 0 || !bcf_hdr_idinfo_exists(hdr, line_type, id)) return BCF_VL_FIXED;
    return bcf_hdr_id2length(hdr, line_type, id);
}

/**
 * @brief Reports whether a length class ties a field's values to the record's ALT list.
 * @param[in] len_class Length class from length_class()
 * @return True if the field's Number is A, R, or G
 */
static bool alt_indexed(int len_class) {
    return len_class == BCF_VL_A || len_class == BCF_VL_R || len_class == BCF_VL_G;
}

/**
 * @brief Splits a record's FORMAT and sample columns into the keys and values worth preserving.
 *
 * A key the writer already emits is dropped rather than duplicated. An ALT-indexed value is kept
 * verbatim and subset to the emitted allele at write time, since the ALT ordinal it is subset
 * against belongs to the variant rather than to the record. The sample may stop short of the key
 * list, which the VCF spec allows, so a value it never supplied is written as missing.
 * @param[in] fmt_col FORMAT column of one record
 * @param[in] sample_col Sample column of one record
 * @param[out] keys Preserved keys, each prefixed with ':'
 * @param[out] vals Preserved values parallel to keys, each prefixed with ':'
 */
static void keep_format(const std::string & fmt_col, const std::string & sample_col,
        std::string & keys, std::string & vals) {
    if (fmt_col == ".") return; // record declares no FORMAT keys at all
    const std::vector<std::string> src_keys = split(fmt_col, ':');
    const std::vector<std::string> src_vals = split(sample_col, ':');
    for (size_t i = 0; i < src_keys.size(); i++) {
        const std::string & key = src_keys[i];
        if (key.empty() || FIXED_FMT_KEY_SET.count(key)) continue;
        keys += ":" + key;
        vals += ":" + (i < src_vals.size() && !src_vals[i].empty() ? src_vals[i] : ".");
    }
}

/**
 * @brief Retains the header lines declaring every field the summary VCF writer preserves.
 *
 * Every declaration is propagated verbatim, an ALT-indexed one included: the record it describes
 * here is biallelic, so Number=A already means one value, Number=R two, and Number=G that record's
 * own genotype count, which is what subsetting leaves behind. The length class is recorded
 * alongside, since the writer subsets without a header to look it up in. PASS and the fixed FORMAT
 * keys are declared by the writer itself.
 * @param[in] hdr Header of the VCF being parsed
 * @param[out] src Store receiving one hdr_keys/hdr_lines entry per retained declaration, and one
 *             info_lens/fmt_lens entry per ALT-indexed field
 */
static void retain_header_lines(const bcf_hdr_t * hdr, std::shared_ptr<srcRecords> src) {
    kstring_t line = {0, 0, NULL};
    for (int i = 0; i < hdr->nhrec; i++) {
        const bcf_hrec_t * hrec = hdr->hrec[i];
        const int line_type = hrec->type;
        if (line_type != BCF_HL_FLT && line_type != BCF_HL_INFO && line_type != BCF_HL_FMT)
            continue;

        std::string id;
        for (int j = 0; j < hrec->nkeys; j++)
            if (std::string(hrec->keys[j]) == "ID") id = hrec->vals[j];
        if (id.empty()) continue;
        if (line_type == BCF_HL_FLT && id == "PASS") continue;
        if (line_type == BCF_HL_FMT && FIXED_FMT_KEY_SET.count(id)) continue;

        // a FILTER declares no Number, so only INFO and FORMAT can be ALT-indexed
        int len_class = line_type == BCF_HL_FLT ? BCF_VL_FIXED : length_class(hdr, line_type, id);
        if (alt_indexed(len_class))
            (line_type == BCF_HL_INFO ? src->info_lens : src->fmt_lens)[id] = len_class;

        line.l = 0;
        if (bcf_hrec_format(hrec, &line) < 0) continue;
        std::string text(line.s, line.l);
        while (!text.empty() && text.back() == '\n') text.pop_back();
        src->hdr_keys.push_back((line_type == BCF_HL_FLT ? "FILTER" :
                line_type == BCF_HL_INFO ? "INFO" : "FORMAT") + ("/" + id));
        src->hdr_lines.push_back(text);
    }
    free(line.s);
}

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

/** @brief Column value stood in for a variant whose source record was not retained. */
static const std::string DOT = ".";
static const std::string PASS = "PASS";
static const std::string NO_EXTRA_FMT = "";

/**
 * @brief Stores one record's preserved columns at its 0-based ordinal within the input VCF.
 *
 * Ordinals arrive in increasing order but not contiguously, since a record dropped before
 * retention still consumes one; the gap it leaves is filled with an empty entry.
 * @param[in] rec_idx 0-based ordinal of the record within its input VCF
 * @param[in] id ID column, verbatim
 * @param[in] qual QUAL column, verbatim
 * @param[in] filter FILTER column, verbatim
 * @param[in] info INFO column, preserved fields only
 * @param[in] fmt_keys Preserved FORMAT keys, each prefixed with ':'
 * @param[in] fmt_vals Sample values parallel to fmt_keys, each prefixed with ':'
 */
void srcRecords::add(int rec_idx, const std::string & id, const std::string & qual,
        const std::string & filter, const std::string & info,
        const std::string & fmt_keys, const std::string & fmt_vals) {
    if (rec_idx >= int(this->ids.size())) {
        this->ids.resize(rec_idx+1);
        this->quals.resize(rec_idx+1);
        this->filters.resize(rec_idx+1);
        this->infos.resize(rec_idx+1);
        this->fmt_keys.resize(rec_idx+1);
        this->fmt_vals.resize(rec_idx+1);
    }
    this->ids[rec_idx] = id;
    this->quals[rec_idx] = qual;
    this->filters[rec_idx] = filter;
    this->infos[rec_idx] = info;
    this->fmt_keys[rec_idx] = fmt_keys;
    this->fmt_vals[rec_idx] = fmt_vals;
}


/**
 * @brief Releases the spare capacity that appending record by record left behind.
 *
 * Growing a vector one record at a time leaves it holding up to twice the memory its contents
 * need. The store is read for the rest of the run and never appended to again once its input is
 * parsed, so that slack is returned before the alignment phase allocates against it.
 */
void srcRecords::shrink() {
    this->ids.shrink_to_fit();
    this->quals.shrink_to_fit();
    this->filters.shrink_to_fit();
    this->infos.shrink_to_fit();
    this->fmt_keys.shrink_to_fit();
    this->fmt_vals.shrink_to_fit();
}


/**
 * @brief Returns a variant's entry in a retained column, or the fallback if it has none.
 *
 * A retained column is empty at every ordinal whose record was dropped before retention, and the
 * whole store is absent for a container built without parsing a VCF, so both stand in as "no
 * source record". An empty FORMAT key list is indistinguishable from an absent one, which is
 * correct: both render as no appended keys.
 * @param[in] column Retained column to read, or nullptr when nothing was retained
 * @param[in] vi Variant index in this container
 * @param[in] fallback Value to return when the variant has no retained entry
 * @return The retained entry, or fallback
 */
const std::string & ctgVariants::src_field(const std::vector<std::string> * column, int vi,
        const std::string & fallback) const {
    if (column == nullptr || vi < 0 || vi >= this->n) return fallback;
    int rec_idx = this->rec_idxs[vi];
    if (rec_idx < 0 || rec_idx >= int(column->size()) || (*column)[rec_idx].empty())
        return fallback;
    return (*column)[rec_idx];
}

/**
 * @brief Returns the source record's ID column, or "." if none was retained.
 * @param[in] vi Variant index in this container
 * @return The ID column
 */
const std::string & ctgVariants::src_id(int vi) const {
    return this->src_field(this->src_recs ? &this->src_recs->ids : nullptr, vi, DOT);
}

/**
 * @brief Returns the source record's QUAL column, or "." if none was retained.
 * @param[in] vi Variant index in this container
 * @return The QUAL column
 */
const std::string & ctgVariants::src_qual(int vi) const {
    return this->src_field(this->src_recs ? &this->src_recs->quals : nullptr, vi, DOT);
}

/**
 * @brief Returns the source record's FILTER column, or "PASS" if none was retained.
 * @param[in] vi Variant index in this container
 * @return The FILTER column
 */
const std::string & ctgVariants::src_filter(int vi) const {
    return this->src_field(this->src_recs ? &this->src_recs->filters : nullptr, vi, PASS);
}

/**
 * @brief Subsets one ALT-indexed value list to the alleles this variant's output record carries.
 *
 * An output record is biallelic, carrying the reference and the single ALT the variant normalized
 * from, so a Number=A list collapses to that ALT's element, a Number=R list to the reference's and
 * that ALT's, and a Number=G list to the genotypes over those two alleles. Genotype g(j,k) with
 * j <= k sits at index k(k+1)/2 + j, which gives the three diploid entries; a haploid genotype
 * list holds one entry per allele instead.
 * @param[in] len_class Length class of the field (BCF_VL_A, BCF_VL_R, or BCF_VL_G)
 * @param[in] alt_idx 1-based ALT ordinal the variant derives from (-1 = unknown)
 * @param[in] ploidy Variant ploidy, which sets the shape of a Number=G list
 * @param[in] value Comma-separated value list, as the source record wrote it
 * @param[out] subset The kept elements, comma-separated; untouched when the subset is undefined
 * @return True if the subset could be taken, false if the ordinal or the list length forbids it
 */
static bool subset_value(int len_class, int alt_idx, ploidy_t ploidy, const std::string & value,
        std::string & subset) {
    if (alt_idx < 1) return false; // no ALT ordinal to index the list with
    std::vector<int> keep;
    switch (len_class) {
        case BCF_VL_A: keep = {alt_idx - 1}; break;
        case BCF_VL_R: keep = {0, alt_idx}; break;
        case BCF_VL_G: keep = ploidy == PLOIDY_HAPLOID ? std::vector<int>{0, alt_idx} :
                std::vector<int>{0, alt_idx*(alt_idx+1)/2, alt_idx*(alt_idx+1)/2 + alt_idx};
                break;
        default: return false;
    }
    const std::vector<std::string> vals = split(value, ',');
    std::string kept;
    for (int i : keep) {
        if (i >= int(vals.size())) return false; // list is shorter than the ALT list it indexes
        kept += kept.empty() ? vals[i] : "," + vals[i];
    }
    subset = kept;
    return true;
}

/**
 * @brief Returns the source record's INFO column subset to this variant's allele, or ".".
 *
 * A field whose values index the original ALT list is subset to the one ALT this variant emits;
 * one that cannot be subset is left out, since INFO expresses a missing field by its absence.
 * @param[in] vi Variant index in this container
 * @return The INFO column, or "." if no source record was retained or no field survived
 */
std::string ctgVariants::src_info(int vi) const {
    const std::string & info =
            this->src_field(this->src_recs ? &this->src_recs->infos : nullptr, vi, DOT);
    if (this->src_recs == nullptr || info == DOT || this->src_recs->info_lens.empty()) return info;

    std::string kept;
    for (const std::string & field : split(info, ';')) {
        if (field.empty()) continue;
        const size_t eq = field.find('=');
        const auto len = this->src_recs->info_lens.find(field.substr(0, eq));
        std::string subset = field;
        if (len != this->src_recs->info_lens.end()) {
            if (eq == std::string::npos) continue; // ALT-indexed but valueless: nothing to subset
            if (!subset_value(len->second, this->alt_idxs[vi], this->ploidies[vi],
                    field.substr(eq+1), subset)) continue;
            subset = field.substr(0, eq) + "=" + subset;
        }
        kept += kept.empty() ? subset : ";" + subset;
    }
    return kept.empty() ? DOT : kept;
}

/**
 * @brief Returns the source record's preserved FORMAT keys, each prefixed with ':'.
 * @param[in] vi Variant index in this container
 * @return The keys to append to the fixed FORMAT list, or "" if there are none
 */
const std::string & ctgVariants::src_fmt_keys(int vi) const {
    return this->src_field(this->src_recs ? &this->src_recs->fmt_keys : nullptr, vi, NO_EXTRA_FMT);
}

/**
 * @brief Returns the source sample's FORMAT values subset to this variant's allele.
 *
 * A field whose values index the original ALT list is subset to the one ALT this variant emits;
 * one that cannot be subset is reported as missing, so the value list stays parallel to the key
 * list src_fmt_keys() returns.
 * @param[in] vi Variant index in this container
 * @return The values to append to the fixed sample fields, each prefixed with ':', or "" if none
 */
std::string ctgVariants::src_fmt_vals(int vi) const {
    const std::string & vals =
            this->src_field(this->src_recs ? &this->src_recs->fmt_vals : nullptr, vi, NO_EXTRA_FMT);
    if (this->src_recs == nullptr || vals.empty() || this->src_recs->fmt_lens.empty()) return vals;

    // both columns are ':'-prefixed, so splitting each leaves a leading empty field to skip
    const std::vector<std::string> keys = split(this->src_fmt_keys(vi), ':');
    const std::vector<std::string> src_vals = split(vals, ':');
    std::string kept;
    for (size_t i = 1; i < keys.size(); i++) {
        const std::string & val = i < src_vals.size() ? src_vals[i] : DOT;
        const auto len = this->src_recs->fmt_lens.find(keys[i]);
        std::string subset;
        if (len == this->src_recs->fmt_lens.end()) subset = val;
        else if (!subset_value(len->second, this->alt_idxs[vi], this->ploidies[vi], val, subset))
            subset = DOT; // no element can be tied to the emitted allele, so none is reported
        kept += ":" + subset;
    }
    return kept;
}


/* Writing back the retained source columns *******************************************************/

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

/**
 * @brief Returns a retained value split on ',', with each VCF missing value typed as htslib's.
 * @param[in] value One field's value text, a comma-separated list for a multi-valued field
 * @return One element per listed value
 */
template <typename T>
static std::vector<T> typed_values(const std::string & value) {
    std::vector<T> values;
    for (const std::string & item : split(value, ',')) {
        T typed;
        if (item == ".") set_missing(typed);
        else if (std::is_same<T, float>::value) typed = T(std::stod(item));
        else typed = T(std::stol(item));
        values.push_back(typed);
    }
    return values;
}

/**
 * @brief Sets one INFO field of a record from its retained text, typed by the output header.
 *
 * The declaration this re-types against is the input's own, copied into the output header, so a
 * field always round-trips through the type it was declared with.
 * @param[in] hdr Summary VCF header, which must declare the field
 * @param[in,out] rec Record to set the field on
 * @param[in] field One ';'-separated INFO entry, either "KEY=VALUE" or a bare Flag key
 * @throws ERROR The field is not declared in the header, or htslib rejects its value
 */
static void set_info_field(const bcf_hdr_t* hdr, bcf1_t* rec, const std::string & field) {
    const size_t eq = field.find('=');
    const std::string key = field.substr(0, eq);
    const std::string value = eq == std::string::npos ? "" : field.substr(eq+1);

    const int id = bcf_hdr_id2int(hdr, BCF_DT_ID, key.data());
    if (id < 0 || !bcf_hdr_idinfo_exists(hdr, BCF_HL_INFO, id))
        ERROR("INFO/%s is not declared in the summary VCF header", key.data());

    int status = 0;
    switch (bcf_hdr_id2type(hdr, BCF_HL_INFO, id)) {
        case BCF_HT_FLAG:
            status = bcf_update_info_flag(hdr, rec, key.data(), NULL, 1);
            break;
        case BCF_HT_INT: {
            const std::vector<int32_t> ints = typed_values<int32_t>(value);
            status = bcf_update_info_int32(hdr, rec, key.data(), ints.data(), int(ints.size()));
            break;
        }
        case BCF_HT_REAL: {
            const std::vector<float> reals = typed_values<float>(value);
            status = bcf_update_info_float(hdr, rec, key.data(), reals.data(), int(reals.size()));
            break;
        }
        default: // BCF_HT_STR, which htslib stores as the value text verbatim
            status = bcf_update_info_string(hdr, rec, key.data(), value.data());
            break;
    }
    if (status < 0) ERROR("Failed to set INFO/%s on summary VCF record", key.data());
}

/**
 * @brief Sets a record's preserved INFO fields, leaving INFO empty when none were retained.
 * @param[in] hdr Summary VCF header, which must declare every preserved field
 * @param[in,out] rec Record to set the fields on
 * @param[in] info Retained INFO column, or "." when the record carried none
 * @throws ERROR A field is not declared in the header, or htslib rejects its value
 */
static void set_record_info(const bcf_hdr_t* hdr, bcf1_t* rec, const std::string & info) {
    if (info == ".") return;
    for (const std::string & field : split(info, ';')) {
        if (!field.empty()) set_info_field(hdr, rec, field);
    }
}

/**
 * @brief Sets a record's QUAL from its retained text, reporting missing for an unscored record.
 * @param[in,out] rec Record to set QUAL on
 * @param[in] qual Retained QUAL column, or "." when the record reported no quality
 */
static void set_record_qual(bcf1_t* rec, const std::string & qual) {
    if (qual == ".") bcf_float_set_missing(rec->qual);
    else rec->qual = float(std::stod(qual));
}

/**
 * @brief Sets a record's ID from its retained text, leaving it missing for an unnamed record.
 * @param[in] hdr Summary VCF header
 * @param[in,out] rec Record to set the ID on
 * @param[in] id Retained ID column, or "." when the record carried no ID
 * @param[in] ctg Contig name, for the error message
 * @param[in] pos Variant position, for the error message
 * @throws ERROR htslib rejects the ID
 */
static void set_record_id(const bcf_hdr_t* hdr, bcf1_t* rec, const std::string & id,
        const std::string & ctg, int pos) {
    if (id == ".") return; // bcf_update_id() with no value already leaves ID missing
    if (bcf_update_id(hdr, rec, id.data()) < 0)
        ERROR("Failed to set ID on summary VCF record at %s:%d", ctg.data(), pos);
}

/**
 * @brief Sets a record's FILTER from its retained text, preserving a non-PASS filter verbatim.
 *
 * A GA4GH consumer reads a non-PASS FILTER on an evaluated record as a filtered call and demotes
 * it, but rewriting the caller's own filter to PASS would misreport the input, so it is carried
 * over as-is and the consequence is documented on write_summary_vcf().
 * @param[in] hdr Summary VCF header, which must declare every filter named
 * @param[in,out] rec Record to set FILTER on
 * @param[in] filter Retained FILTER column, ';'-separated, or "." when the record was unfiltered
 * @param[in] ctg Contig name, for the error message
 * @param[in] pos Variant position, for the error message
 * @throws ERROR A named filter is not declared in the header, or htslib rejects it
 */
static void set_record_filters(const bcf_hdr_t* hdr, bcf1_t* rec, const std::string & filter,
        const std::string & ctg, int pos) {
    if (filter == ".") return; // an unfiltered record declares no filter at all
    for (const std::string & name : split(filter, ';')) {
        if (name.empty()) continue;
        const int id = bcf_hdr_id2int(hdr, BCF_DT_ID, name.data());
        if (id < 0 || !bcf_hdr_idinfo_exists(hdr, BCF_HL_FLT, id))
            ERROR("FILTER '%s' is not declared in the summary VCF header", name.data());
        if (bcf_add_filter(hdr, rec, id) < 0)
            ERROR("Failed to set FILTER on summary VCF record at %s:%d", ctg.data(), pos);
    }
}

/**
 * @brief Sets the FORMAT fields carried over from the record's owning sample.
 *
 * One FORMAT key list serves both samples, so the sample that does not own the record has no
 * values for these keys and reports each as missing.
 * @param[in] hdr Summary VCF header, which must declare every preserved field
 * @param[in,out] rec Record whose fixed FORMAT fields are already set
 * @param[in] owner The owning sample's values, whose src_keys and src_vals are read
 * @param[in] query_owns True if the owner is the query sample, which is written second
 * @throws ERROR A field is not declared in the header, or htslib rejects its values
 */
static void set_source_formats(const bcf_hdr_t* hdr, bcf1_t* rec, const sample_fields & owner,
        bool query_owns) {
    if (owner.src_keys.empty()) return;
    const std::vector<std::string> keys = split(owner.src_keys, ':');
    const std::vector<std::string> vals = split(owner.src_vals, ':');

    // both lists lead with an empty element, since every entry is written ':'-prefixed
    for (size_t i = 1; i < keys.size(); i++) {
        const std::string & key = keys[i];
        const std::string value = i < vals.size() ? vals[i] : ".";
        const int id = bcf_hdr_id2int(hdr, BCF_DT_ID, key.data());
        if (id < 0 || !bcf_hdr_idinfo_exists(hdr, BCF_HL_FMT, id))
            ERROR("FORMAT/%s is not declared in the summary VCF header", key.data());

        int status = 0;
        switch (bcf_hdr_id2type(hdr, BCF_HL_FMT, id)) {
            case BCF_HT_INT: {
                const std::vector<int32_t> mine = typed_values<int32_t>(value);
                const std::vector<int32_t> padded = query_owns ?
                        pad_per_allele(std::vector<int32_t>{}, mine) :
                        pad_per_allele(mine, std::vector<int32_t>{});
                status = bcf_update_format_int32(hdr, rec, key.data(), padded.data(),
                        int(padded.size()));
                break;
            }
            case BCF_HT_REAL: {
                const std::vector<float> mine = typed_values<float>(value);
                const std::vector<float> padded = query_owns ?
                        pad_per_allele(std::vector<float>{}, mine) :
                        pad_per_allele(mine, std::vector<float>{});
                status = bcf_update_format_float(hdr, rec, key.data(), padded.data(),
                        int(padded.size()));
                break;
            }
            default: { // BCF_HT_STR, which htslib stores as one string per sample
                const char* values[2] = {query_owns ? "." : value.data(),
                                         query_owns ? value.data() : "."};
                status = bcf_update_format_string(hdr, rec, key.data(), values, 2);
                break;
            }
        }
        if (status < 0) ERROR("Failed to set FORMAT/%s on summary VCF record", key.data());
    }
}

/**************************************************************************************************/

/**
 * @brief Builds the summary VCF header, declaring every FORMAT field and the TRUTH/QUERY samples.
 *
 * Every FILTER, INFO, and FORMAT field carried over from an input is declared too, so that the
 * output describes itself and htslib accepts the values written against it. A record's columns
 * come from whichever callset owns it, so both callsets contribute declarations; where the two
 * disagree about an ID the query's wins, as it does for the columns of a matched record.
 * @param[in] contigs Contig names, in the order records are written
 * @param[in] lengths Contig lengths, parallel to contigs
 * @param[in] src_recs Retained source records of each callset, any of which may be nullptr
 * @return Header owning its own memory, to be released by the caller with bcf_hdr_destroy()
 * @throws ERROR The header cannot be allocated, a header line htslib rejects, a sample htslib
 *         rejects, or a header htslib cannot synchronize
 */
bcf_hdr_t* summary_vcf_header(const std::vector<std::string> & contigs,
        const std::vector<int> & lengths,
        const EnumArray<callset_t, std::shared_ptr<srcRecords>, CALLSET_SLOTS> & src_recs /* = {} */) {

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

    // A source record that parsing split -- a complex variant into an INS and a DEL, or a het-alt
    // into one entry per ALT -- becomes several records here, and each resolves back to the same
    // source record, so each carries an identical copy of its preserved columns. That is intended,
    // but a reader cannot infer it from the file, and summing a count-like preserved field over
    // these records double-counts the one source value.
    lines.push_back("##vcfdistPreservedFields=<Description=\"ID, QUAL, FILTER, INFO, and the "
            "non-fixed FORMAT fields are carried over from whichever callset owns each record. "
            "One source record may yield several records here (a complex variant is split into an "
            "INS and a DEL, a het-alt into one record per ALT), each repeating the same preserved "
            "values, so summing a count-like field over records double-counts the source value. "
            "Number=A/R/G values index the source ALT list, so each is subset to the one allele "
            "its record carries.\">");

    // declare the fields carried over from the inputs, PASS excluded since it is declared above
    std::unordered_set<std::string> declared = {"FILTER/PASS"};
    for (callset_t c : EnumRange<callset_t, CALLSET_SLOTS>{}) {
        if (src_recs[c] == nullptr) continue;
        const srcRecords & src = *src_recs[c];
        for (size_t i = 0; i < src.hdr_lines.size(); i++) {
            if (declared.insert(src.hdr_keys[i]).second) lines.push_back(src.hdr_lines[i]);
        }
    }

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
 * @brief Sets the fixed VCF fields (CHROM, POS, ID, REF, ALT, QUAL, FILTER, INFO) of one record.
 *
 * ID, QUAL, FILTER, and INFO come from the source record, as POS, REF, and ALT already do. INFO
 * fields whose values index the source ALT list are subset to the one allele this record emits. A
 * variant with no retained source record falls back to the placeholders '.', '.', PASS, and '.'
 * that this writer emitted before any field was preserved.
 * @param[in] hdr Summary VCF header, which must declare this contig and every preserved field
 * @param[in,out] rec Cleared record to fill
 * @param[in] ref Reference FASTA data for retrieving flanking bases for indels
 * @param[in] ctg Contig name
 * @param[in] idx Variant index in this container
 * @throws ERROR The contig is not declared in the header
 * @throws ERROR An INS/DEL sits at the contig start (0-based pos 0), leaving no preceding base to anchor
 * @throws ERROR The variant type is not TYPE_SUB, TYPE_INS, or TYPE_DEL
 * @throws ERROR htslib rejects the record's ID, FILTER, alleles, or any preserved INFO field
 */
void ctgVariants::set_var_record(const bcf_hdr_t* hdr, bcf1_t* rec,
        std::shared_ptr<fastaData> ref, const std::string & ctg, int idx) const {

    rec->rid = bcf_hdr_name2id(hdr, ctg.data());
    if (rec->rid < 0) ERROR("Contig '%s' is not declared in the summary VCF header", ctg.data());
    set_record_qual(rec, this->src_qual(idx));
    set_record_id(hdr, rec, this->src_id(idx), ctg, this->poss[idx]);
    set_record_filters(hdr, rec, this->src_filter(idx), ctg, this->poss[idx]);

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

    set_record_info(hdr, rec, this->src_info(idx));
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
 * @param[in] owns_record If true, this sample's source record supplies the record's preserved
 *        columns, so its own FORMAT fields are carried over; the other sample reports them missing
 * @return This sample's FORMAT values, htslib-encoded
 */
sample_fields ctgVariants::var_sample_fields(int vi, int sc_idx, int phase_block,
        bool phase_switch, bool phase_flip, bool query /* = false */,
        bool owns_record /* = false */) const {

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

    // the record's FORMAT key list comes from its owner alone, so the other sample adds none
    if (owns_record) {
        fields.src_keys = this->src_fmt_keys(vi);
        fields.src_vals = this->src_fmt_vals(vi);
    }
    return fields;
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

    // only the sample owning the record carries preserved fields, so at most one side has any
    set_source_formats(hdr, rec, truth.src_keys.empty() ? query : truth, truth.src_keys.empty());
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

    // source column retention, reusing one render buffer across every record
    kstring_t rec_str = {0, 0, NULL};

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

    // the summary VCF is the only consumer of the retained source columns, and it is written on
    // every run, so retention is unconditional; the declarations are read once, from this header
    variant_data->src_recs = std::shared_ptr<srcRecords>(new srcRecords());
    retain_header_lines(hdr, variant_data->src_recs);

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
        variant_data->variants[HAP1][ctgnames[i]]->src_recs = variant_data->src_recs;
        variant_data->variants[HAP2][ctgnames[i]]->src_recs = variant_data->src_recs;
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

        // retain the columns the summary VCF writer cannot derive from the internal arrays;
        // rendering the whole record reuses htslib's own formatting instead of re-deriving the
        // text of each field, and indexing by record ordinal shares one copy across the entries
        // a multi-allelic or complex record splits into
        if (variant_data->src_recs != nullptr) {
            rec_str.l = 0;
            if (vcf_format(hdr, rec, &rec_str) < 0)
                ERROR("Failed to format %s VCF record at %s:%lld",
                        callset_strs[callset].data(), ctg.data(), (long long)rec->pos);
            std::string rendered(rec_str.s, rec_str.l);
            while (!rendered.empty() && rendered.back() == '\n') rendered.pop_back();
            std::vector<std::string> cols = split(rendered, '\t');
            if (int(cols.size()) < VCF_COLS)
                ERROR("Expected %d columns but found %d in %s VCF record at %s:%lld",
                        int(VCF_COLS), int(cols.size()), callset_strs[callset].data(),
                        ctg.data(), (long long)rec->pos);
            std::string fmt_keys, fmt_vals;
            keep_format(cols[FORMAT_COL], cols[SAMPLE_COL], fmt_keys, fmt_vals);
            variant_data->src_recs->add(n-1, cols[ID_COL], cols[QUAL_COL], cols[FILTER_COL],
                    cols[INFO_COL], fmt_keys, fmt_vals);
        }

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

    if (variant_data->src_recs != nullptr) variant_data->src_recs->shrink();

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
    free(rec_str.s);
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
