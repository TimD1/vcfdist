/**
 * @file bed.cpp
 * @brief BED file loading, interval storage, and contig intersection utilities.
 */
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
 * The first three columns are read and the remainder are ignored.
 * @param[in] bed_fn The BED filename.
 * @throws ERROR if the BED file cannot be opened.
 * @throws ERROR if a start or stop coordinate is empty, not numeric in full, or too large.
 */
bedData::bedData(const std::string & bed_fn) {

    // fail if file doesn't exist
    auto bed_fp = fopen(bed_fn.data(), "r");
    if (bed_fp == NULL) {
        ERROR("Failed to open BED file '%s'", bed_fn.data());
    }
    fclose(bed_fp);

    std::ifstream bed(bed_fn);
    std::string region;
    int line = 0;
    while (getline(bed, region)) {
        line++;
        std::stringstream ss(region);
        std::string contig, start, stop;
        getline(ss, contig, '\t');
        getline(ss, start, '\t');
        getline(ss, stop, '\t');
        // parsed into locals so that the reported field does not depend on evaluation order
        const int start_pos = parse_coord(start, bed_fn, line);
        const int stop_pos = parse_coord(stop, bed_fn, line);
        this->add(contig, start_pos, stop_pos);
    }
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
 * @brief Returns BED location type for a variant interval.
 *
 * @param[in] contig The contig containing the variant.
 * @param[in] start The 0-based inclusive start position of the variant.
 * @param[in] stop The 0-based exclusive end position of the variant.
 * @param[in] type The type of the variant.
 * @return One of: BED_INSIDE, BED_OUTSIDE, BED_BORDER, BED_OFFCTG.
 * @throws ERROR if the variant stop precedes the variant start.
 */
int bedData::contains(std::string contig, const int & start, const int & stop, const int & type) {

    if (!g.bed_exists) return BED_INSIDE;

    if (stop < start)
        ERROR("Invalid region %s:%d-%d in BED contains", contig.data(), start, stop);

    // contig not in BED
    if (this->regions.find(contig) == this->regions.end()) return BED_OFFCTG;

    // variant before/after all BED regions
    if (stop <= this->regions[contig].starts[0]) return BED_OUTSIDE;
    if (start >= this->regions[contig].stops[ 
            this->regions[contig].stops.size()-1]) return BED_OUTSIDE;

    // get indices of variant within bed regions list
    int start_idx = std::upper_bound(
            this->regions[contig].starts.begin(),
            this->regions[contig].starts.end(),
            start) - this->regions[contig].starts.begin() - 1;
    int stop_idx = std::lower_bound(
            this->regions[contig].stops.begin(),
            this->regions[contig].stops.end(),
            stop) - this->regions[contig].stops.begin();

    // variant must be partially in region, other index off end
    if (start_idx < 0) return BED_BORDER;
    if (stop_idx >= int(this->regions[contig].stops.size())) return BED_BORDER;

    // variant in middle
    if (stop_idx == start_idx) {
        // don't allow INS exactly at region end
        if (type == TYPE_INS && start == this->regions[contig].stops[stop_idx]-1)
            return BED_BORDER;
        return BED_INSIDE;
    }
    if (stop_idx == start_idx + 1) {
        int next_region_start = this->regions[contig].starts[stop_idx];
        int prev_region_stop = this->regions[contig].stops[start_idx];
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
 * @brief Intersects reference FASTA, query VCF, truth VCF, and optional BED regions, retaining only common contigs.
 *
 * @param[in] query_ptr A pointer to the query variantData.
 * @param[in] truth_ptr A pointer to the truth variantData.
 * @param[in] ref_ptr A pointer to the reference fastaData.
 * @throws WARNING if contigs in either VCF are not present in either the other VCF or BED file.
 * @throws WARNING if corresponding contigs in the truth and query VCFs observed differing ploidies.
 * @throws ERROR if a contig to be evaluated is not present in the reference FASTA.
 */
void intersect_contigs(
        std::shared_ptr<variantData> query_ptr,
        std::shared_ptr<variantData> truth_ptr,
        std::shared_ptr<fastaData> ref_ptr) {
    if (g.verbosity >= 1) INFO(" ");
    if (g.verbosity >= 1) INFO("  Checking contigs:");

    if (g.bed_exists) { // use BED to determine contigs

        // remove all extraneous contigs in query VCF not in BED
        std::vector<std::string>::iterator itr = query_ptr->contigs.begin();
        while (itr != query_ptr->contigs.end()) { // query
            if (std::find(g.bed.contigs.begin(), g.bed.contigs.end(),
                        *itr) == g.bed.contigs.end()) {
                query_ptr->lengths.erase(query_ptr->lengths.begin() + 
                        (itr - query_ptr->contigs.begin()));
                query_ptr->observed_ploidies.erase(query_ptr->observed_ploidies.begin() +
                        (itr - query_ptr->contigs.begin()));
                query_ptr->variants[HAP1].erase(*itr);
                query_ptr->variants[HAP2].erase(*itr);
                std::string dropped_ctg = *itr; // save name, erase() invalidates itr
                itr = query_ptr->contigs.erase(itr);
                if (g.verbosity >= 2)
                    WARN("Ignoring %s from QUERY VCF, not in BED file.", dropped_ctg.data());
            } else ++itr;
        }
        // remove all extraneous contigs in truth VCF not in BED
        itr = truth_ptr->contigs.begin();
        while (itr != truth_ptr->contigs.end()) { // truth
            if (std::find(g.bed.contigs.begin(), g.bed.contigs.end(),
                        *itr) == g.bed.contigs.end()) {
                truth_ptr->lengths.erase(truth_ptr->lengths.begin() + 
                        (itr - truth_ptr->contigs.begin()));
                truth_ptr->observed_ploidies.erase(truth_ptr->observed_ploidies.begin() +
                        (itr - truth_ptr->contigs.begin()));
                truth_ptr->variants[HAP1].erase(*itr);
                truth_ptr->variants[HAP2].erase(*itr);
                std::string dropped_ctg = *itr; // save name, erase() invalidates itr
                itr = truth_ptr->contigs.erase(itr);
                if (g.verbosity >= 2)
                    WARN("Ignoring %s from TRUTH VCF, not in BED file.", dropped_ctg.data());
            } else ++itr;
        }
        // remove all extraneous contigs in ref FASTA not in BED
        auto itr2 = ref_ptr->fasta.begin();
        while (itr2 != ref_ptr->fasta.end()) { // fasta
            if (std::find(g.bed.contigs.begin(), g.bed.contigs.end(),
                        itr2->first) == g.bed.contigs.end()) {
                itr2 = ref_ptr->fasta.erase(itr2);
            } else itr2++;
        }

        // warn if list of truth and query contigs are not the same
        for (std::string ctg : query_ptr->contigs) {
            if (std::find(truth_ptr->contigs.begin(), 
                        truth_ptr->contigs.end(), ctg) == truth_ptr->contigs.end())
                WARN("Contig '%s' found in query VCF but not truth VCF."
                     " All query variants on '%s' will be false positives.", ctg.data(), ctg.data());
        }
        for (std::string ctg : truth_ptr->contigs) {
            if (std::find(query_ptr->contigs.begin(), 
                        query_ptr->contigs.end(), ctg) == query_ptr->contigs.end())
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

        // ensure fasta contains all contigs
        for (std::string ctg : truth_ptr->contigs) {
            if (ref_ptr->fasta.find(ctg) == ref_ptr->fasta.end())
                ERROR("Contig '%s' found in truth VCF but not reference FASTA. Please provide BED file.", ctg.data());
        }

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

        // remove extra contigs from ref FASTA
        auto itr = ref_ptr->fasta.begin();
        while (itr != ref_ptr->fasta.end()) {
            if (std::find(truth_ptr->contigs.begin(), truth_ptr->contigs.end(),
                        itr->first) == truth_ptr->contigs.end()) {
                itr = ref_ptr->fasta.erase(itr);
            } else itr++;
        }
    }

    // verify the observed ploidies match for all truth/query contigs
    for (int i = 0; i < int(truth_ptr->contigs.size()); i++) {
        std::string ctg = truth_ptr->contigs[i];
        int query_ctg_idx = std::find(query_ptr->contigs.begin(),
                query_ptr->contigs.end(), ctg) - query_ptr->contigs.begin();
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
