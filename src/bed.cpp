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
    if (regions.find(contig) == regions.end()) {
        regions[contig] = contigRegions();
    }
    regions[contig].starts.push_back(start);
    regions[contig].stops.push_back(stop);
}

void bedData::check() {
}

int bedData::contains(std::string contig, const int & start, const int & stop) {

    printf("@%s %d-%d\n", contig.data(), start, stop);
    if (stop < start)
        ERROR("Invalid region %s:%d-%d in BED contains", contig.data(), start, stop);

    // contig not in BED
    if (this->regions.find(contig) == this->regions.end()) {
        printf("no contig\n");
        return OUTSIDE;
    }

    // variant before/after all BED regions
    if (stop <= this->regions[contig].starts[0]) {
        printf("stop before first region\n");
        return OUTSIDE;
    }
    if (start >= this->regions[contig].stops[ 
            this->regions[contig].stops.size()-1]) {
        printf("start after last region\n");
        return OUTSIDE;
    }

    // get indices of variant within bed regions list
    int start_idx = std::upper_bound(
            this->regions[contig].starts.begin(),
            this->regions[contig].starts.end(),
            start) - this->regions[contig].starts.begin() - 1;
    int stop_idx = std::lower_bound(
            this->regions[contig].stops.begin(),
            this->regions[contig].stops.end(),
            stop) - this->regions[contig].stops.begin();

    // we know stop > starts[0], and start < starts[0]
    if (start_idx < 0) {
        printf("start before first region, stop must be after first region starts.\n");
        return BORDER;
    }
    if (stop_idx >= this->regions[contig].stops.size()) {
        printf("stop after last region, start must be before last region starts.\n");
        return BORDER;
    }

    int region_start = this->regions[contig].starts[start_idx];
    int region_stop = this->regions[contig].stops[stop_idx];
    printf("start idx: %d   stop_idx: %d\n", start_idx, stop_idx);
    printf("region_start: %d   region_stop: %d\n", region_start, region_stop);
    if (stop_idx == start_idx)
        return INSIDE;
    if (stop_idx == start_idx + 1) {
        int next_region_start = this->regions[contig].starts[stop_idx];
        int prev_region_stop = this->regions[contig].stops[start_idx];
        if (start >= prev_region_stop && stop <= next_region_start) 
            return OUTSIDE;
        return BORDER;
    }
    return BORDER;
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
