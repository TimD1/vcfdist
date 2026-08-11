/**
 * @file bed.cpp
 * @brief BED file loading, interval storage, and contig intersection utilities.
 */
#include <filesystem>
#include <utility>

#include "htslib/hts.h"
#include "htslib/kstring.h"

#include "bed.h"
#include "print.h"


/* File-local helpers *****************************************************************************/

/**
 * @brief Parses one BED coordinate field, naming the offending line if it is not a valid int.
 *
 * The whole field must be numeric, so a partially-numeric coordinate such as "8bp" is rejected
 * rather than silently truncated to the digits that precede the suffix.
 * @param[in] coord The coordinate field, exactly as read from the BED file.
 * @param[in] bed_fn The BED filename, reported on failure.
 * @param[in] line The 1-based line the field was read from, reported on failure.
 * @return The parsed coordinate.
 * @throws ERROR if the field is empty, not numeric in full, or too large for an int.
 */
static int parse_coord(const std::string & coord, const std::string & bed_fn, const int & line) {
    int pos = 0;
    bool valid = false;
    try {
        size_t len = 0;
        pos = std::stoi(coord, &len);
        valid = len == coord.size();
    } catch (const std::exception &) { // non-numeric or out-of-range, reported below
    }
    if (!valid) {
        ERROR("Invalid coordinate '%s' on line %d of BED file '%s'",
                coord.data(), line, bed_fn.data());
    }
    return pos;
}


/* bedData ****************************************************************************************/

/**
 * @brief Constructs a bedData object by reading intervals from a BED file.
 *
 * The first three columns are read and the remainder are ignored. Reading goes through htslib, so
 * plain, gzip, and bgzip encodings are all accepted; the encoding is detected from the file's
 * leading bytes rather than from its extension.
 *
 * The loaded regions are optionally normalized, then validated with check() either way. Rejecting
 * an unsorted or overlapping file outright is right for an evaluation region whose malformation
 * silently changes every denominator; normalizing first is right for a third-party region set we
 * neither author nor control. Normalizing does not weaken the validation, since disorder and
 * overlap are the only malformations it repairs -- a flipped or zero-length interval still errors.
 * @param[in] bed_fn The BED filename.
 * @param[in] normalize Whether to sort and merge the regions before validating them.
 * @throws ERROR if the BED file cannot be opened.
 * @throws ERROR if a line cannot be read, which a truncated or corrupt compressed file looks like.
 * @throws ERROR if a start or stop coordinate is empty, not numeric in full, or too large.
 * @throws ERROR if the regions fail check().
 * @throws WARNING if the regions are normalized and were not already sorted.
 * @throws WARNING if adjacent regions share a boundary, which normalizing would have merged.
 */
bedData::bedData(const std::string & bed_fn, bool normalize) {

    this->filename = bed_fn;

    // fail if file doesn't exist, or is compressed in a way htslib cannot decode
    htsFile* bed_fp = hts_open(bed_fn.data(), "r");
    if (bed_fp == NULL) {
        ERROR("Failed to open BED file '%s'", bed_fn.data());
    }

    kstring_t region = KS_INITIALIZE;
    int line = 0;
    int len = 0;
    while ((len = hts_getline(bed_fp, '\n', &region)) >= 0) {
        line++;
        // hts_getline strips the terminator, so the line is the record and nothing else
        std::stringstream ss(std::string(region.s == NULL ? "" : region.s, region.l));
        std::string contig, start, stop;
        getline(ss, contig, '\t');
        getline(ss, start, '\t');
        getline(ss, stop, '\t');
        // parsed into locals so that the reported field does not depend on evaluation order
        const int start_pos = parse_coord(start, bed_fn, line);
        const int stop_pos = parse_coord(stop, bed_fn, line);
        this->add(contig, start_pos, stop_pos);
    }
    // -1 is end-of-file; anything lower is a read failure, which must not look like a short file
    if (len < -1) {
        ERROR("Failed to read line %d of BED file '%s'", line+1, bed_fn.data());
    }
    ks_free(&region);
    hts_close(bed_fp);

    // normalizing repairs disorder and overlap; check() still has to reject what it cannot
    if (normalize) this->normalize();
    this->check();
}

/**
 * @brief Adds a single [start, stop) region on a named contig to this bedData.
 *
 * @param[in] contig The name of the contig.
 * @param[in] start The 0-based inclusive start position of the interval.
 * @param[in] stop The 0-based exclusive end position of the interval.
 */
void bedData::add(const std::string & contig, const int & start, const int & stop) {
    // add a new contig if needed
    if (this->regions.find(contig) == this->regions.end()) {
        this->regions[contig] = contigRegions();
        this->contigs.push_back(contig);
    }
    // add the region
    this->regions[contig].starts.push_back(start);
    this->regions[contig].stops.push_back(stop);
    this->regions[contig].n++;
    this->size += stop-start;
}

/**
 * @brief Validates that all BED intervals are sorted and non-overlapping.
 *
 * @throws ERROR if a BED region is flipped (stop < start), zero-length, unsorted, or overlapping.
 * @throws WARNING if adjacent BED regions share a boundary and should be merged.
 */
void bedData::check() {
    for (size_t ctg_idx = 0; ctg_idx < this->contigs.size(); ctg_idx++) {
        int prev_start = -1;
        int prev_stop = -1;
        for (int region_idx = 0; region_idx < this->regions[this->contigs[ctg_idx]].n; region_idx++) {

            // check this region is positive size
            int start = this->regions[this->contigs[ctg_idx]].starts[region_idx];
            int stop = this->regions[this->contigs[ctg_idx]].stops[region_idx];
            if (stop < start) ERROR("BED region %s:%d-%d stop precedes start.", 
                    this->contigs[ctg_idx].data(), start, stop);
            if (stop == start) ERROR("BED region %s:%d-%d length zero.",
                    this->contigs[ctg_idx].data(), start, stop);

            // check for overlaps/etc
            if (region_idx) { // not first region on contig

                if (stop < prev_start)
                    ERROR("BED is unsorted; region %s:%d-%d precedes %s:%d-%d.",
                        this->contigs[ctg_idx].data(), prev_start, prev_stop, 
                        this->contigs[ctg_idx].data(), start, stop);
                if (start < prev_stop)
                    ERROR("BED overlap detected: regions %s:%d-%d and %s:%d-%d.",
                        this->contigs[ctg_idx].data(), prev_start, prev_stop, 
                        this->contigs[ctg_idx].data(), start, stop);
                if (prev_stop == start) 
                    WARN("BED regions %s:%d-%d and %s:%d-%d should be merged.",
                        this->contigs[ctg_idx].data(), prev_start, prev_stop, 
                        this->contigs[ctg_idx].data(), start, stop);
            }

            prev_start = start;
            prev_stop = stop;
        }
    }
}

/**
 * @brief Sorts each contig's intervals by start and merges those that overlap or abut.
 *
 * contains() locates a variant with two binary searches over the start and stop lists, so an
 * unsorted or overlapping region set makes it return wrong answers rather than merely being
 * untidy. The total size is recomputed from the merged intervals, since bases covered by more than
 * one input interval were counted once per interval by add().
 *
 * A flipped or zero-length interval is left exactly as it was found, since neither is repairable
 * by sorting or merging; check() rejects both.
 * @throws WARNING if the intervals on any contig were not already sorted.
 */
void bedData::normalize() {

    long merged_size = 0;
    int coalesced = 0;
    bool was_unsorted = false;
    for (const std::string & contig : this->contigs) {
        contigRegions & ctg_regions = this->regions[contig];

        // sort by start, then by stop, so that a region nested in another follows it
        std::vector< std::pair<int, int> > intervals;
        for (int i = 0; i < ctg_regions.n; i++)
            intervals.push_back(std::make_pair(ctg_regions.starts[i], ctg_regions.stops[i]));
        // reported rather than silently corrected, since an unsorted region set means the file is
        // not what check() would have accepted
        if (!std::is_sorted(intervals.begin(), intervals.end())) was_unsorted = true;
        std::sort(intervals.begin(), intervals.end());

        std::vector<int> starts, stops;
        for (const auto & [start, stop] : intervals) {
            // starts ascend, so any interval reaching the running stop extends it instead
            if (!starts.empty() && start <= stops.back()) {
                stops.back() = std::max(stops.back(), stop);
            } else {
                starts.push_back(start);
                stops.push_back(stop);
            }
        }
        for (size_t i = 0; i < starts.size(); i++) merged_size += stops[i] - starts[i];

        coalesced += ctg_regions.n - int(starts.size());
        ctg_regions.starts = std::move(starts);
        ctg_regions.stops = std::move(stops);
        ctg_regions.n = int(ctg_regions.starts.size());
    }
    this->size = merged_size;

    if (was_unsorted)
        WARN("Regions in BED file '%s' were unsorted, and have been sorted.",
                this->filename.data());
    if (g.verbosity >= 2 && coalesced)
        INFO("Merged %d overlapping or adjacent regions in BED file '%s'.",
                coalesced, this->filename.data());
}

/**
 * @brief Returns BED location type for a variant interval.
 *
 * @param[in] contig The contig containing the variant.
 * @param[in] start The 0-based inclusive start position of the variant.
 * @param[in] stop The 0-based exclusive end position of the variant.
 * @param[in] type The type of the variant.
 * @return One of: BED_INSIDE, BED_OUTSIDE, BED_BORDER, BED_OFFCTG.
 * @throws ERROR if the variant stop precedes the variant start.
 */
bedloc_t bedData::contains(std::string contig, const int & start, const int & stop,
        const edittype_t & type) {

    if (stop < start)
        ERROR("Invalid region %s:%d-%d in BED contains", contig.data(), start, stop);

    // contig not in BED
    if (this->regions.find(contig) == this->regions.end()) return BED_OFFCTG;

    // get indices of variant within bed regions list
    int start_idx = std::upper_bound(
            this->regions[contig].starts.begin(),
            this->regions[contig].starts.end(),
            start) - this->regions[contig].starts.begin() - 1;
    int stop_idx = std::lower_bound(
            this->regions[contig].stops.begin(),
            this->regions[contig].stops.end(),
            stop) - this->regions[contig].stops.begin();

    return this->classify(contig, start, stop, type, start_idx, stop_idx);
}

/**
 * @brief Returns BED location type for a variant already located within a contig's intervals.
 *
 * Locating a variant and classifying it are kept apart so that a caller locating it some other way
 * -- a cursor advanced across position-sorted variants, say -- shares this decision tree instead of
 * reimplementing it. No search is performed here; both indices are supplied by the caller.
 * @param[in] contig The contig containing the variant, which must be present in this bedData.
 * @param[in] start The 0-based inclusive start position of the variant.
 * @param[in] stop The 0-based exclusive end position of the variant.
 * @param[in] type The type of the variant.
 * @param[in] start_idx Index of the last region starting at or before start, or -1 if there is none.
 * @param[in] stop_idx Index of the first region stopping at or after stop, or the region count.
 * @return One of: BED_INSIDE, BED_OUTSIDE, BED_BORDER.
 */
bedloc_t bedData::classify(const std::string & contig, const int & start, const int & stop,
        const edittype_t & type, const int & start_idx, const int & stop_idx) {

    const contigRegions & ctg_regions = this->regions.at(contig);

    // variant before/after all BED regions; these precede the index tests below because a variant
    // left of the first region has start_idx -1, which would otherwise read as BED_BORDER
    if (stop <= ctg_regions.starts[0]) return BED_OUTSIDE;
    if (start >= ctg_regions.stops.back()) return BED_OUTSIDE;

    // variant must be partially in region, other index off end
    if (start_idx < 0) return BED_BORDER;
    if (stop_idx >= int(ctg_regions.stops.size())) return BED_BORDER;

    // variant in middle
    if (stop_idx == start_idx) {
        // don't allow INS exactly at region end
        if (type == TYPE_INS && start == ctg_regions.stops[stop_idx]-1)
            return BED_BORDER;
        return BED_INSIDE;
    }
    if (stop_idx == start_idx + 1) {
        int next_region_start = ctg_regions.starts[stop_idx];
        int prev_region_stop = ctg_regions.stops[start_idx];
        if (start >= prev_region_stop && stop <= next_region_start)
            return BED_OUTSIDE; // between
        return BED_BORDER; // both
    }
    return BED_BORDER; // spans multiple regions
}

/**
 * @brief Returns a formatted string listing all stored BED regions by contig.
 */
bedData::operator std::string() const {
    std::string bed_regions = "";
    for (const auto &[contig, region_list]: this->regions) {
        bed_regions += contig + ":\n";
        for (int i = 0; i < region_list.n; i++) {
            bed_regions += "\t" + std::to_string(region_list.starts[i]) + "-" + \
                           std::to_string(region_list.stops[i]) + "\n";
        }
    }
    return bed_regions;
}


/* Stratification manifest ************************************************************************/

/**
 * @brief Resolves one manifest BED path against the directory holding the manifest.
 *
 * The GIAB stratification manifests name their region sets with paths relative to themselves, so
 * resolving against the working directory would require every manifest to be rewritten before use.
 * @param[in] path The BED path exactly as written in the manifest.
 * @param[in] strat_tsv_fn The manifest filename, whose parent directory relative paths resolve to.
 * @return The path to open, unchanged if it was already absolute.
 */
static std::string resolve_strat_path(const std::string & path, const std::string & strat_tsv_fn) {
    if (!path.empty() && path[0] == '/') return path;
    return parent_path(strat_tsv_fn) + path;
}

/**
 * @brief Loads every region set named by a stratification manifest TSV, in manifest order.
 *
 * The manifest is hap.py-compatible: two tab-separated columns naming a stratum and its BED file,
 * with any further columns ignored. Blank lines and lines beginning with '#' are skipped. Reading
 * goes through htslib, so the manifest and each region set may be plain, gzip, or bgzip encoded.
 *
 * Each region set is normalized at load, so it is sorted and merged before being validated; these
 * are third-party files we neither author nor control, and contains() would return wrong answers on
 * an unsorted or overlapping one. Normalizing does not weaken the validation: a malformation
 * sorting and merging cannot repair still errors. Both outputs are cleared first, so a repeated
 * call replaces the strata rather than appending to them.
 * @param[in] strat_tsv_fn The stratification manifest filename.
 * @param[out] strat_names Stratum names, in manifest order.
 * @param[out] strata Parsed stratum regions, parallel to strat_names.
 * @throws ERROR if the manifest cannot be opened.
 * @throws ERROR if a manifest line cannot be read, which a truncated compressed file looks like.
 * @throws ERROR if a line names fewer than two fields.
 * @throws ERROR if a stratum name repeats, which would make an output row ambiguous.
 * @throws ERROR if a stratum is named '*', which is reserved for the all-regions row.
 * @throws ERROR if a stratum's BED file cannot be opened.
 * @throws WARNING if the manifest names no region sets at all.
 */
void load_strata(const std::string & strat_tsv_fn, std::vector<std::string> & strat_names,
        std::vector<bedData> & strata) {

    strat_names.clear();
    strata.clear();

    htsFile* strat_fp = hts_open(strat_tsv_fn.data(), "r");
    if (strat_fp == NULL) {
        ERROR("Failed to open stratification manifest '%s'", strat_tsv_fn.data());
    }

    kstring_t manifest_line = KS_INITIALIZE;
    int line = 0;
    int len = 0;
    while ((len = hts_getline(strat_fp, '\n', &manifest_line)) >= 0) {
        line++;
        // hts_getline strips the terminator, so the line is the record and nothing else
        const std::string text(manifest_line.s == NULL ? "" : manifest_line.s, manifest_line.l);
        if (text.empty() || text[0] == '#') continue;

        std::stringstream ss(text);
        std::string name, path;
        getline(ss, name, '\t');
        getline(ss, path, '\t');
        if (name.empty() || path.empty()) {
            ERROR("Line %d of stratification manifest '%s' does not name both a stratum and a"
                    " BED file", line, strat_tsv_fn.data());
        }
        if (name == "*") {
            ERROR("Stratum name '*' on line %d of stratification manifest '%s' is reserved for the"
                    " all-regions row", line, strat_tsv_fn.data());
        }
        if (std::find(strat_names.begin(), strat_names.end(), name) != strat_names.end()) {
            ERROR("Duplicate stratum name '%s' on line %d of stratification manifest '%s'",
                    name.data(), line, strat_tsv_fn.data());
        }

        // an unopenable region set names both paths, since a relative path that resolved somewhere
        // unintended is indistinguishable from a missing file when only one of the two is reported
        const std::string bed_fn = resolve_strat_path(path, strat_tsv_fn);
        htsFile* bed_fp = hts_open(bed_fn.data(), "r");
        if (bed_fp == NULL) {
            const std::string abs_fn = std::filesystem::absolute(bed_fn).string();
            ERROR("Failed to open BED file '%s' for stratum '%s' on line %d of stratification"
                    " manifest '%s', resolved to '%s'",
                    path.data(), name.data(), line, strat_tsv_fn.data(), abs_fn.data());
        }
        hts_close(bed_fp);

        strat_names.push_back(name);
        strata.push_back(bedData(bed_fn, true));
    }
    // -1 is end-of-file; anything lower is a read failure, which must not look like a short file
    if (len < -1) {
        ERROR("Failed to read line %d of stratification manifest '%s'",
                line+1, strat_tsv_fn.data());
    }
    ks_free(&manifest_line);
    hts_close(strat_fp);

    if (strat_names.empty()) {
        WARN("Stratification manifest '%s' names no region sets", strat_tsv_fn.data());
    } else if (g.verbosity >= 1) {
        INFO("Loaded %d stratification region sets from '%s'",
                int(strat_names.size()), strat_tsv_fn.data());
    }
}

/**
 * @brief Warns for stratification region sets sharing no contig with the reference FASTA.
 *
 * A region set naming contigs the reference does not have contributes to no output row but a row of
 * zeroes, which reads as a genuine result. The realistic cause is an assembly or contig-naming
 * mismatch -- 'chr20'-prefixed region sets against '20'-named reference contigs, or GRCh37 sets
 * against a GRCh38 run -- so every set failing at once is diagnostic rather than incidental and is
 * reported as one message naming that cause. A set that shares a contig but contains no evaluated
 * variant is not warned about; that is a legitimate result.
 * @param[in] ref_ptr A pointer to the reference fastaData.
 * @throws WARNING naming each region set that shares no contig with the reference.
 * @throws WARNING naming the likely cause instead, if no region set shares a contig at all.
 */
void check_strata_contigs(const std::shared_ptr<fastaData> & ref_ptr) {

    std::vector<std::string> unmatched;
    for (size_t i = 0; i < g.strata.size(); i++) {
        bool shares_contig = false;
        for (const std::string & ctg : g.strata[i].contigs) {
            if (ref_ptr->fasta.find(ctg) != ref_ptr->fasta.end()) {
                shares_contig = true;
                break;
            }
        }
        if (!shares_contig) unmatched.push_back(g.strat_names[i]);
    }
    if (unmatched.empty()) return;

    if (unmatched.size() == g.strata.size()) {
        WARN("No stratification region set in '%s' shares a contig with reference FASTA '%s';"
                " check that both use the same reference assembly and contig naming"
                " (such as 'chr20' rather than '20')",
                g.strat_tsv_fn.data(), g.ref_fasta_fn.data());
    } else {
        for (const std::string & name : unmatched) {
            WARN("Stratification region set '%s' shares no contig with reference FASTA '%s'",
                    name.data(), g.ref_fasta_fn.data());
        }
    }
}


/* BED helper functions ***************************************************************************/

/**
 * @brief Renders a set of observed ploidies as a brace-delimited, comma-separated list.
 * @param[in] ploidies Ploidies observed on one contig in one callset
 * @return The ploidies in ascending order, e.g. "{1}" or "{1,2}"
 */
static std::string ploidy_set_str(const std::set<int> & ploidies) {
    std::string str = "{";
    for (auto itr = ploidies.begin(); itr != ploidies.end(); itr++) {
        if (itr != ploidies.begin()) str += ",";
        str += std::to_string(*itr);
    }
    return str + "}";
}

/**
 * @brief Reconciles the contigs of the reference FASTA, query VCF, truth VCF, and optional BED.
 *
 * Every contig either VCF carries is retained, along with its reference sequence, so that variants
 * on it remain reportable; contigs outside the BED are simply never evaluated, since
 * parse_variants() has already discarded their variants as BED_OFFCTG. Contigs the BED names but a
 * VCF lacks are injected empty, so that both callsets cover every evaluated region.
 *
 * @param[in,out] query_ptr A pointer to the query variantData, gaining any missing BED contigs.
 * @param[in,out] truth_ptr A pointer to the truth variantData, gaining any missing BED contigs.
 * @param[in] ref_ptr A pointer to the reference fastaData.
 * @throws WARNING if a BED contig is present in one VCF but not the other.
 * @throws WARNING if corresponding evaluated contigs in the truth and query VCFs observed differing
 *         ploidies.
 * @throws ERROR if a contig in either VCF is not present in the reference FASTA.
 * @throws ERROR if a contig in the BED file is not present in the reference FASTA.
 */
void intersect_contigs(
        std::shared_ptr<variantData> query_ptr,
        std::shared_ptr<variantData> truth_ptr,
        std::shared_ptr<fastaData> ref_ptr) {
    if (g.verbosity >= 1) INFO(" ");
    if (g.verbosity >= 1) INFO("  Checking contigs:");

    // Every contig in either callset is retained, so every one of them needs reference sequence:
    // set_var_record() anchors an INS/DEL in ref->fasta when writing a variant back out.
    const std::vector< std::pair<callset_t, std::shared_ptr<variantData> > > callsets =
            {{QUERY, query_ptr}, {TRUTH, truth_ptr}};
    for (const auto & [callset, vcf_ptr] : callsets) {
        for (const std::string & ctg : vcf_ptr->contigs) {
            if (ref_ptr->fasta.find(ctg) == ref_ptr->fasta.end())
                ERROR("Contig '%s' found in %s VCF but not reference FASTA.",
                        ctg.data(), callset_strs[callset].data());
        }
    }

    if (g.bed_exists) { // use BED to determine contigs

        // Warn if the truth and query contig lists differ, over the BED contigs alone: nothing is
        // evaluated outside them, so a one-sided contig there yields neither FPs nor FNs to warn
        // about. Runs before the injection below, which would otherwise mask every difference.
        for (const std::string & ctg : g.bed.contigs) {
            bool in_query = std::find(query_ptr->contigs.begin(),
                    query_ptr->contigs.end(), ctg) != query_ptr->contigs.end();
            bool in_truth = std::find(truth_ptr->contigs.begin(),
                    truth_ptr->contigs.end(), ctg) != truth_ptr->contigs.end();
            if (in_query && !in_truth)
                WARN("Contig '%s' found in query VCF but not truth VCF."
                     " All query variants on '%s' will be false positives.", ctg.data(), ctg.data());
            if (in_truth && !in_query)
                WARN("Contig '%s' found in truth VCF but not query VCF."
                     " All truth variants on '%s' will be false negatives.", ctg.data(), ctg.data());
        }

        // ensure all inputs contain required contigs (even if empty)
        for (std::string ctg : g.bed.contigs) {
            if (ref_ptr->fasta.find(ctg) == ref_ptr->fasta.end())
                ERROR("Contig '%s' found in BED but not reference FASTA.", ctg.data());
            if (std::find(query_ptr->contigs.begin(), 
                        query_ptr->contigs.end(), ctg) == query_ptr->contigs.end()) {
                INFO("Contig '%s' found in BED but not query VCF.", ctg.data());
                query_ptr->variants[HAP1][ctg] = 
                        std::shared_ptr<ctgVariants>(new ctgVariants(ctg));
                query_ptr->variants[HAP2][ctg] = 
                        std::shared_ptr<ctgVariants>(new ctgVariants(ctg));
                query_ptr->contigs.push_back(ctg);
                query_ptr->lengths.push_back(ref_ptr->lengths.at(ctg));
                query_ptr->observed_ploidies.push_back({});
            }
            if (std::find(truth_ptr->contigs.begin(),
                        truth_ptr->contigs.end(), ctg) == truth_ptr->contigs.end()) {
                INFO("Contig '%s' found in BED but not truth VCF.", ctg.data());
                truth_ptr->variants[HAP1][ctg] = 
                        std::shared_ptr<ctgVariants>(new ctgVariants(ctg));
                truth_ptr->variants[HAP2][ctg] = 
                        std::shared_ptr<ctgVariants>(new ctgVariants(ctg));
                truth_ptr->contigs.push_back(ctg);
                truth_ptr->lengths.push_back(ref_ptr->lengths.at(ctg));
                truth_ptr->observed_ploidies.push_back({});
            }
        }

    } else { // use truth VCF to determine contigs

        // ensure query/truth VCFs contain the same contigs (even if devoid of variants)
        for (int i = 0; i < int(query_ptr->contigs.size()); i++) {
            std::string ctg = query_ptr->contigs[i];
            if (std::find(truth_ptr->contigs.begin(), 
                        truth_ptr->contigs.end(), ctg) == truth_ptr->contigs.end()) {
                WARN("Contig '%s' found in query VCF but not truth VCF."
                     " All query variants on '%s' will be false positives.", ctg.data(), ctg.data());
                truth_ptr->variants[HAP1][ctg] = 
                        std::shared_ptr<ctgVariants>(new ctgVariants(ctg));
                truth_ptr->variants[HAP2][ctg] = 
                        std::shared_ptr<ctgVariants>(new ctgVariants(ctg));
                truth_ptr->contigs.push_back(ctg);
                truth_ptr->lengths.push_back(ref_ptr->lengths.at(ctg));
                truth_ptr->observed_ploidies.push_back({});
            }
        }
        for (int i = 0; i < int(truth_ptr->contigs.size()); i++) {
            std::string ctg = truth_ptr->contigs[i];
            if (std::find(query_ptr->contigs.begin(), 
                        query_ptr->contigs.end(), ctg) == query_ptr->contigs.end()) {
                WARN("Contig '%s' found in truth VCF but not query VCF."
                     " All truth variants on '%s' will be false negatives.", ctg.data(), ctg.data());
                query_ptr->variants[HAP1][ctg] = 
                        std::shared_ptr<ctgVariants>(new ctgVariants(ctg));
                query_ptr->variants[HAP2][ctg] = 
                        std::shared_ptr<ctgVariants>(new ctgVariants(ctg));
                query_ptr->contigs.push_back(ctg);
                query_ptr->lengths.push_back(ref_ptr->lengths.at(ctg));
                query_ptr->observed_ploidies.push_back({});
            }
        }
    }

    // verify the observed ploidies match for all evaluated truth/query contigs
    for (int i = 0; i < int(truth_ptr->contigs.size()); i++) {
        std::string ctg = truth_ptr->contigs[i];

        // Ploidy is recorded before parse_variants() applies the BED filter, so a contig outside
        // the BED carries ploidies no evaluation consults; disagreement there is not actionable.
        if (g.bed_exists && std::find(g.bed.contigs.begin(), g.bed.contigs.end(), ctg)
                == g.bed.contigs.end()) continue;

        // Unreachable: the skip above drops the BED-absent contigs, and the injection above pairs
        // every remaining one. Guarded anyway, since indexing on a failed find() would read out of
        // bounds rather than report anything.
        auto query_itr = std::find(query_ptr->contigs.begin(), query_ptr->contigs.end(), ctg);
        if (query_itr == query_ptr->contigs.end()) continue;
        int query_ctg_idx = query_itr - query_ptr->contigs.begin();
        int truth_ctg_idx = i;
        const std::set<int> & truth_ploidies = truth_ptr->observed_ploidies[truth_ctg_idx];
        const std::set<int> & query_ploidies = query_ptr->observed_ploidies[query_ctg_idx];

        // a contig injected empty observed no ploidy at all, which is not a disagreement
        if (truth_ploidies.empty() || query_ploidies.empty()) continue;
        if (truth_ploidies != query_ploidies) {
            WARN("%s contig '%s' has ploidies %s and %s contig '%s' has ploidies %s",
                    callset_strs[TRUTH].data(), ctg.data(), ploidy_set_str(truth_ploidies).data(),
                    callset_strs[QUERY].data(), ctg.data(), ploidy_set_str(query_ploidies).data());
        }
    }

    if (g.verbosity >= 1) INFO("    All contig checks passed!");
}
