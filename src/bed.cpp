#include "bed.h"
#include "print.h"

bedData::bedData(const std::string & bed_fn) {

    // fail if file doesn't exist
    auto bed_fp = fopen(bed_fn.data(), "r");
    if (bed_fp == NULL) {
        ERROR("Failed to open BED file '%s'", bed_fn.data());
    }
    fclose(bed_fp);

    std::ifstream bed(bed_fn);
    std::string region;
    while (getline(bed, region)) {
        std::stringstream ss(region);
        std::string contig, start, stop;
        getline(ss, contig, '\t');
        getline(ss, start, '\t');
        getline(ss, stop, '\t');
        this->add(contig, std::stoi(start), std::stoi(stop));
    }
}

void bedData::add(const std::string & contig, const int & start, const int & stop) {
    if (this->regions.find(contig) == this->regions.end()) {
        this->regions[contig] = contigRegions();
        this->contigs.push_back(contig);
    }
    this->regions[contig].starts.push_back(start);
    this->regions[contig].stops.push_back(stop);
    this->regions[contig].n++;
    this->size += stop-start;
}

void bedData::check() {
    for (size_t i = 0; i < this->contigs.size(); i++) {
        int prev_start = -1;
        int prev_stop = -1;
        for (int j = 0; j < this->regions[this->contigs[i]].n; j++) {

            // check this region is positive size
            int start = this->regions[this->contigs[i]].starts[j];
            int stop = this->regions[this->contigs[i]].stops[j];
            if (stop < start) ERROR("BED region %s:%d-%d stop precedes start.", 
                    this->contigs[i].data(), start, stop);
            if (stop == start) ERROR("BED region %s:%d-%d length zero.",
                    this->contigs[i].data(), start, stop);

            // check for overlaps/etc
            if (j) { // not first region on contig

                if (stop < prev_start)
                    ERROR("BED is unsorted; region %s:%d-%d precedes %s:%d-%d.",
                        this->contigs[i].data(), prev_start, prev_stop, 
                        this->contigs[i].data(), start, stop);
                if (start < prev_stop)
                    ERROR("BED overlap detected: regions %s:%d-%d and %s:%d-%d.",
                        this->contigs[i].data(), prev_start, prev_stop, 
                        this->contigs[i].data(), start, stop);
                if (prev_stop == start) 
                    WARN("BED regions %s:%d-%d and %s:%d-%d should be merged.",
                        this->contigs[i].data(), prev_start, prev_stop, 
                        this->contigs[i].data(), start, stop);
            }

            prev_start = start;
            prev_stop = stop;
        }
    }
}

int bedData::contains(std::string contig, const int & start, const int & stop) {

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
    if (stop_idx == start_idx)
        return BED_INSIDE;
    if (stop_idx == start_idx + 1) {
        int next_region_start = this->regions[contig].starts[stop_idx];
        int prev_region_stop = this->regions[contig].stops[start_idx];
        if (start >= prev_region_stop && stop <= next_region_start) 
            return BED_OUTSIDE; // between
        return BED_BORDER; // both
    }
    return BED_BORDER; // spans multiple regions
}

bedData::operator std::string() const {
    std::string bed_regions = "";
    for (const auto &[contig, region_list]: this->regions) {
        bed_regions += contig + ":\n";
        for (size_t i = 0; i < region_list.starts.size(); i++) {
            bed_regions += "\t" + std::to_string(region_list.starts[i]) + "-" + \
                           std::to_string(region_list.stops[i]) + "\n";
        }
    }
    return bed_regions;
}


/******************************************************************************/


void check_contigs(
        std::unique_ptr<variantData> & query_ptr,
        std::unique_ptr<variantData> & truth_ptr,
        std::shared_ptr<fastaData> ref_ptr) {
    if (g.verbosity >= 1) INFO("Checking contigs");

    if (g.bed_exists) { // use BED to determine contigs

        // remove all extraneous contigs not in BED
        std::vector<std::string>::iterator itr = query_ptr->contigs.begin();
        while (itr != query_ptr->contigs.end()) { // query
            if (std::find(g.bed.contigs.begin(), g.bed.contigs.end(),
                        *itr) == g.bed.contigs.end()) {
                query_ptr->lengths.erase(query_ptr->lengths.begin() + 
                        (itr - query_ptr->contigs.begin()));
                query_ptr->ctg_variants[HAP1].erase(*itr);
                query_ptr->ctg_variants[HAP2].erase(*itr);
                itr = query_ptr->contigs.erase(itr);
            } else itr++;
        }
        itr = truth_ptr->contigs.begin();
        while (itr != truth_ptr->contigs.end()) { // truth
            if (std::find(g.bed.contigs.begin(), g.bed.contigs.end(),
                        *itr) == g.bed.contigs.end()) {
                truth_ptr->lengths.erase(truth_ptr->lengths.begin() + 
                        (itr - truth_ptr->contigs.begin()));
                truth_ptr->ctg_variants[HAP1].erase(*itr);
                truth_ptr->ctg_variants[HAP2].erase(*itr);
                itr = truth_ptr->contigs.erase(itr);
            } else itr++;
        }
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
                WARN("Contig '%s' found in query VCF but not truth VCF.", ctg.data());
        }
        for (std::string ctg : truth_ptr->contigs) {
            if (std::find(query_ptr->contigs.begin(), 
                        query_ptr->contigs.end(), ctg) == query_ptr->contigs.end())
                WARN("Contig '%s' found in truth VCF but not query VCF.", ctg.data());
        }

        // ensure all inputs contain required contigs (even if empty)
        for (std::string ctg : g.bed.contigs) {
            if (ref_ptr->fasta.find(ctg) == ref_ptr->fasta.end())
                ERROR("Contig '%s' found in BED but not reference FASTA.", ctg.data());
            if (std::find(query_ptr->contigs.begin(), 
                        query_ptr->contigs.end(), ctg) == query_ptr->contigs.end()) {
                INFO("Contig '%s' found in BED but not query VCF.", ctg.data());
                query_ptr->ctg_variants[HAP1][ctg] = 
                        std::shared_ptr<ctgVariants>(new ctgVariants());
                query_ptr->ctg_variants[HAP2][ctg] = 
                        std::shared_ptr<ctgVariants>(new ctgVariants());
                query_ptr->contigs.push_back(ctg);
                query_ptr->lengths.push_back(ref_ptr->lengths.at(ctg));
            }
            if (std::find(truth_ptr->contigs.begin(), 
                        truth_ptr->contigs.end(), ctg) == truth_ptr->contigs.end()) {
                INFO("Contig '%s' found in BED but not truth VCF.", ctg.data());
                truth_ptr->ctg_variants[HAP1][ctg] = 
                        std::shared_ptr<ctgVariants>(new ctgVariants());
                truth_ptr->ctg_variants[HAP2][ctg] = 
                        std::shared_ptr<ctgVariants>(new ctgVariants());
                truth_ptr->contigs.push_back(ctg);
                truth_ptr->lengths.push_back(ref_ptr->lengths.at(ctg));
            }
        }

    } else { // use truth VCF to determine contigs

        // ensure fasta contains all contigs
        for (std::string ctg : truth_ptr->contigs) {
            if (ref_ptr->fasta.find(ctg) == ref_ptr->fasta.end())
                ERROR("Contig '%s' found in truth VCF but not reference FASTA. Please provide BED file.", ctg.data());
        }

        // remove all extraneous contigs
        std::vector<std::string>::iterator itr = query_ptr->contigs.begin();
        while (itr != query_ptr->contigs.end()) { // query
            if (std::find(truth_ptr->contigs.begin(), truth_ptr->contigs.end(),
                        *itr) == truth_ptr->contigs.end()) {
                ERROR("Contig '%s' found in query VCF but not truth VCF. Please provide BED file.", (*itr).data());
                query_ptr->lengths.erase(query_ptr->lengths.begin() + 
                        (itr - query_ptr->contigs.begin()));
                query_ptr->ctg_variants[HAP1].erase(*itr);
                query_ptr->ctg_variants[HAP2].erase(*itr);
                itr = query_ptr->contigs.erase(itr);
            } else itr++;
        }
        auto itr2 = ref_ptr->fasta.begin();
        while (itr2 != ref_ptr->fasta.end()) { // fasta
            if (std::find(truth_ptr->contigs.begin(), truth_ptr->contigs.end(),
                        itr2->first) == truth_ptr->contigs.end()) {
                itr2 = ref_ptr->fasta.erase(itr2);
            } else itr2++;
        }

        // ensure all inputs contain required contigs (even if empty)
        for (std::string ctg : truth_ptr->contigs) {
            if (std::find(query_ptr->contigs.begin(), 
                        query_ptr->contigs.end(), ctg) == query_ptr->contigs.end()) {
                WARN("Contig '%s' found in truth VCF but not query VCF.", ctg.data());
                query_ptr->ctg_variants[HAP1][ctg] = 
                        std::shared_ptr<ctgVariants>(new ctgVariants());
                query_ptr->ctg_variants[HAP2][ctg] = 
                        std::shared_ptr<ctgVariants>(new ctgVariants());
                query_ptr->contigs.push_back(ctg);
                query_ptr->lengths.push_back(ref_ptr->lengths.at(ctg));
            }
        }
    }
}
