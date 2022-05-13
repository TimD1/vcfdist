#ifndef _BED_H_
#define _BED_H_

#include <string>
#include <vector>
#include <unordered_map>
#include <fstream>
#include <sstream>

struct contigRegions {
    std::vector<int> starts;
    std::vector<int> stops;
};

class bedData {
public:


    // constructors
    bedData() {;}
    bedData(const std::string & bed_fn) {
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

    void add(const std::string & contig, const int & start, const int & stop) {
        if (regions.find(contig) == regions.end()) {
            regions[contig] = contigRegions();
        }
        regions[contig].starts.push_back(start);
        regions[contig].stops.push_back(stop);
    }

    bool includes(std::string contig, long pos) {
        return false;
    }

    operator std::string() const {
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

private:
    std::unordered_map<std::string, contigRegions> regions;
};

#endif
