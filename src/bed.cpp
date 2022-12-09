#include "bed.h"
#include "print.h"

bedData::bedData(const std::string & bed_fn) {
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

    if (stop < start)
        ERROR("Invalid region %s:%d-%d in BED contains", contig.data(), start, stop);

    // contig not in BED
    if (this->regions.find(contig) == this->regions.end()) return BED_OFFCTG;

    // variant before/after all BED regions
    if (stop <= this->regions[contig].starts[0]) return BED_OUTSIDE;
    if (start >= this->regions[contig].stops[ 
            this->regions[contig].stops.size()-1]) return BED_OUTSIDE;

    // get indices of variant within bed regions list
    size_t start_idx = std::upper_bound(
            this->regions[contig].starts.begin(),
            this->regions[contig].starts.end(),
            start) - this->regions[contig].starts.begin() - 1;
    size_t stop_idx = std::lower_bound(
            this->regions[contig].stops.begin(),
            this->regions[contig].stops.end(),
            stop) - this->regions[contig].stops.begin();

    // variant must be partially in region, other index off end
    if (start_idx < 0) return BED_BORDER;
    if (stop_idx >= this->regions[contig].stops.size()) return BED_BORDER;

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
