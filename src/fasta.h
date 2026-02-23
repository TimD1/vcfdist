/**
 * @file fasta.h
 * @brief FASTA reference sequence storage and indexing by contig name.
 */
#ifndef _FASTA_H_
#define _FASTA_H_

#include <algorithm>
#include <string>
#include <unordered_map>

// zlib is required for kseq
#include "zlib.h"
#include "htslib/kseq.h"
KSEQ_INIT(int, read);

/**
 * @class fastaData
 * @brief Stores all sequences from a FASTA reference file, indexed by contig name.
 */
class fastaData {
public:
    /**
     * @brief Reads and stores all sequences from an open FASTA file, converting to uppercase.
     * @param[in] ref_fasta_fp Open file pointer to FASTA file (closed after reading)
     */
    fastaData(FILE * ref_fasta_fp) {
        kseq_t * seq = kseq_init(fileno(ref_fasta_fp));
        while (kseq_read(seq) >= 0) {
            this->fasta[seq->name.s] = seq->seq.s;
            std::transform(this->fasta[seq->name.s].begin(), this->fasta[seq->name.s].end(),
                    this->fasta[seq->name.s].begin(), ::toupper);
            this->lengths[seq->name.s] = this->fasta.at(seq->name.s).size();
        }
        kseq_destroy(seq);
        fclose(ref_fasta_fp);
    }

    std::unordered_map<std::string,std::string> fasta;   ///< Map from contig name to uppercase sequence string
    std::unordered_map<std::string,int> lengths;         ///< Map from contig name to sequence length
};

#endif
