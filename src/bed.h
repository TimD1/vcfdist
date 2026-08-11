/**
 * @file bed.h
 * @brief BED file loading, interval storage, and contig intersection utilities.
 */
#ifndef _BED_H_
#define _BED_H_

#include <string>
#include <vector>
#include <unordered_map>
#include <sstream>
#include <algorithm>

#include "defs.h"
#include "fasta.h"
#include "variant.h"

/**
 * @struct contigRegions
 * @brief Stores n sorted [start, stop) intervals present on a given contig in a BED file.
 */
struct contigRegions {
    std::vector<int> starts; ///< 0-based start coordinates for intervals
    std::vector<int> stops; ///< 0-based non-inclusive end coordinates for intervals
    int n = 0; ///< number of intervals
};

/**************************************************************************************************/

/**
 * @class bedData
 * @brief Loads and stores interval information contained in a BED file (first 3 columns).
 *
 * This class is used for defining regions of interest to evaluate (masking ground truth low
 * confidence regions), and therefore does not store additional BED columns.
 */
class bedData {
public:

    // constructors
    /** @brief Constructs an empty bedData object with no regions. */
    bedData() {;}
    /** @brief Constructs a bedData object by reading regions from a BED file, optionally normalizing. */
    bedData(const std::string & bed_fn, bool normalize = false);

    // member functions
    /** @brief Adds a single [start, stop) region on a contig to this bedData. */
    void add(const std::string & contig, const int & start, const int & stop);
    /** @brief Validates that all BED intervals are sorted and non-overlapping. */
    void check();
    /** @brief Sorts each contig's intervals and merges those that overlap or abut. */
    void normalize();
    /** @brief Returns BED location type (BED_INSIDE/OUTSIDE/BORDER/OFFCTG) for a variant interval. */
    bedloc_t contains(std::string contig, const int & start, const int & stop,
            const edittype_t & type);
    /** @brief Returns BED location type for a variant already located within a contig's intervals. */
    bedloc_t classify(const std::string & contig, const int & start, const int & stop,
            const edittype_t & type, const int & start_idx, const int & stop_idx);

    /** @brief Returns a string representation of all stored BED regions. */
    operator std::string() const;

    // member variables
    std::string filename; ///< BED file the regions were read from, empty if built in memory
    long size = 0;  ///< total size of all regions in BED file, in bases
    std::vector<std::string> contigs; ///< list of contigs in BED file
    std::unordered_map<std::string, contigRegions> regions; ///< mapping from contigs to struct storing all intervals on a contig
};

/** @brief Loads every region set named by a stratification manifest TSV, in manifest order. */
void load_strata(const std::string & strat_tsv_fn, std::vector<std::string> & strat_names,
        std::vector<bedData> & strata);

/** @brief Warns for stratification region sets sharing no contig with the reference FASTA. */
void check_strata_contigs(const std::shared_ptr<fastaData> & ref_ptr);

/** @brief Intersects query, truth, and reference contigs with BED regions and retains only common contigs. */
void intersect_contigs(
        std::shared_ptr<variantData> query_ptr,
        std::shared_ptr<variantData> truth_ptr,
        std::shared_ptr<fastaData> ref_ptr);

#endif
