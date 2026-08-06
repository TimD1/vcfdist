/**
 * @file phase.cpp
 * @brief Phase block detection, switch/flip error classification, and phasing summary output.
 */
#include <vector>
#include <algorithm>
#include <stdexcept>
#include <string>
#include <unordered_set>

#include "phase.h"
#include "print.h"
#include "globals.h"


/**************************************************************************************************/

/** @brief Retained records of one contig, indexed by callset (nullptr = that callset retained none). */
typedef std::vector< std::shared_ptr<ctgSideline> > ctg_sidelines;

/** @brief Retained records of every contig, indexed by callset then contig. */
typedef EnumArray<callset_t, std::unordered_map< std::string, std::shared_ptr<ctgSideline> >,
        CALLSET_SLOTS> all_sidelines;

/**
 * @brief Returns each callset's retained records for one contig, nullptr where it retained none.
 * @param[in] all Retained records of every contig, by callset
 * @param[in] ctg Contig to look up
 * @return One entry per callset, in callset order
 */
static ctg_sidelines sidelines_of(const all_sidelines & all, const std::string & ctg) {
    ctg_sidelines side(CALLSETS, nullptr);
    for (callset_t c : EnumRange<callset_t, CALLSET_SLOTS>{}) {
        const auto found = all[c].find(ctg);
        if (found != all[c].end()) side[idx(c)] = found->second;
    }
    return side;
}

/**
 * @brief Writes the retained records of one contig that begin before a position, in position order.
 *
 * The query's records come first where the two callsets tie. A tie against an evaluated variant is
 * resolved the other way, by the caller passing that variant's position as the limit: a retained
 * record belongs to no supercluster, so it has no place within the evaluated ordering of a position.
 * @param[in] out_fp Open file pointer to output VCF
 * @param[in] ctg Contig name
 * @param[in] side Retained records of this contig, by callset
 * @param[in,out] side_ptrs Index of each callset's next unwritten record, advanced as records write
 * @param[in] limit Position (0-based) at or past which records are left for a later call
 */
static void write_sidelined(FILE* out_fp, const std::string & ctg, const ctg_sidelines & side,
        std::vector<int> & side_ptrs, int limit) {
    while (true) {
        int next = -1;
        int next_pos = limit;
        for (int c = 0; c < CALLSETS; c++) {
            if (side[c] == nullptr || side_ptrs[c] >= side[c]->n) continue;
            if (side[c]->poss[side_ptrs[c]] < next_pos) {
                next_pos = side[c]->poss[side_ptrs[c]];
                next = c;
            }
        }
        if (next < 0) return;
        side[next]->print_var(out_fp, ctg, side_ptrs[next], callset_t(next));
        side_ptrs[next]++;
    }
}

/**
 * @brief Returns the contigs holding retained records but no evaluated ones, in input header order.
 *
 * A contig the BED file never mentions is dropped from the evaluated contig list before any
 * analysis runs, so the walk over that list never reaches it, yet its records are retained and
 * still belong in the output. One container per such contig is returned, the query's where both
 * callsets retained records on it, since only the contig's name and length are read from it.
 * @param[in] all Retained records of every contig, by callset
 * @param[in] evaluated Contigs the evaluated walk already covers
 * @return One container per uncovered contig holding records, ordered by input header ordinal
 */
static ctg_sidelines retained_only_contigs(const all_sidelines & all,
        const std::vector<std::string> & evaluated) {
    ctg_sidelines found;
    for (callset_t c : EnumRange<callset_t, CALLSET_SLOTS>{}) {
        for (const auto & [ctg, side] : all[c]) {
            if (side->n == 0) continue;
            if (std::find(evaluated.begin(), evaluated.end(), ctg) != evaluated.end()) continue;
            bool seen = false;
            for (const auto & other : found) seen = seen || other->ctg == ctg;
            if (!seen) found.push_back(side);
        }
    }
    // the map they were collected from is unordered, so the sort is what makes the output stable
    std::sort(found.begin(), found.end(), [](const std::shared_ptr<ctgSideline> & a,
            const std::shared_ptr<ctgSideline> & b) {
        return (a->rid != b->rid) ? a->rid < b->rid : a->ctg < b->ctg;
    });
    return found;
}

/**
 * @brief Writes a summary VCF containing all variants annotated with benchmark metrics.
 * @param[in] out_vcf_fn Output VCF filename
 * @note FORMAT fields include: TP/FP/FN decision, credit score, edit distances, phase info, and flip/switch errors
 * @note One record is written per variant, each sample's GT reporting that callset's own claim; the
 *       per-haplotype fields are lists carrying one value per GT allele
 * @note A het-alt (1|2) source record stays two co-located records, since parsing splits it into
 *       two entries whose alleles normalize independently
 * @note Contigs called by only one callset are included; a contig with no query variants has no
 *       phase block to read, so its truth records are written with PB and BS defaulted
 * @note ID, QUAL, FILTER, INFO, and the appended FORMAT keys come from whichever callset owns the
 *       record, which is the query wherever it calls and the truth only on a pure false negative.
 *       A matched record therefore drops the truth record's INFO and FORMAT, and the callset that
 *       does not own a record writes '.' for each appended FORMAT key
 * @note Number=A/R/G values index the source record's ALT list, while these records carry
 *       normalized, split alleles, so each is subset to the one allele its record emits. The
 *       propagated declaration is rewritten to match: Number=A becomes 1 and Number=R becomes 2,
 *       while Number=G stands, already resolving to the genotype count of the record's ploidy
 * @note A variant retained without being evaluated is held outside the evaluated arrays, so it is
 *       merged in by position rather than walked alongside them. It was excluded before any
 *       comparison ran, so it is never matched: one sample carries its call, reporting BD=N and
 *       nothing else, and the other is entirely missing
 * @note A reason decided per allele can retain one haplotype of a record while the other is
 *       evaluated, in which case the record appears twice: once normalized to the allele that was
 *       evaluated, and once whole, since nothing split the retained one
 * @note A contig absent from the BED file is dropped from the evaluated contig list before any
 *       analysis runs, so its retained records are declared and written after every evaluated
 *       contig rather than in the walk over that list
 * @note The input FILTER is preserved verbatim, on evaluated records included. A retained record
 *       adds the VCFDIST_-prefixed tag naming why it was not evaluated, replacing a lone PASS.
 *       A GA4GH consumer reads a non-PASS FILTER on an evaluated record as a filtered call and
 *       demotes it, turning filtered TPs into FNs and filtered FPs into Ns, so a caller's own
 *       non-PASS filter accepted via --filter will demote those calls downstream. Rewriting it to
 *       PASS would misreport the input, so it is reported here rather than designed around
 * @throws ERROR if the output summary VCF file cannot be opened for writing
 * @throws ERROR if neither callset is selected next while variants remain
 */
void phaseblockData::write_summary_vcf(std::string out_vcf_fn) {

    // VCF header
    if (g.verbosity >= 1) INFO("  Writing summary VCF to '%s'", out_vcf_fn.data());
    FILE* out_vcf = fopen(out_vcf_fn.data(), "w");
    if (out_vcf == NULL) {
        ERROR("Failed to open summary VCF file '%s'", out_vcf_fn.data());
    }
    const std::chrono::time_point<std::chrono::system_clock> now{std::chrono::system_clock::now()};
    time_t tt = std::chrono::system_clock::to_time_t(now);
    tm local_time = *localtime(&tt);
    fprintf(out_vcf, "##fileformat=VCFv4.2\n");
    fprintf(out_vcf, "##fileDate=%04d%02d%02d\n", local_time.tm_year + 1900, 
            local_time.tm_mon + 1, local_time.tm_mday);
    fprintf(out_vcf, "##CL=%s\n", g.cmd.data());
    for (size_t i = 0; i < this->contigs.size(); i++) {
        fprintf(out_vcf, "##contig=<ID=%s,length=%d>\n",
                this->contigs[i].data(), this->lengths[i]);
    }

    // A contig absent from the BED file is dropped before evaluation, but its records are retained
    // and written below, so it is declared here: a record on an undeclared contig is not valid VCF.
    // These follow the evaluated contigs in the header because their records follow in the body.
    const ctg_sidelines retained_only =
            retained_only_contigs(this->callset_sidelined, this->contigs);
    for (const std::shared_ptr<ctgSideline> & side : retained_only) {
        fprintf(out_vcf, "##contig=<ID=%s,length=%d>\n", side->ctg.data(), side->length);
    }
    fprintf(out_vcf, "##FILTER=<ID=PASS,Description=\"All filters passed\">\n");

    // Declare the tag each retained record carries. The IDs are prefixed so that no input FILTER
    // ID can collide with one, and uppercase to match the convention for FILTER IDs.
    std::unordered_set<std::string> declared = {"FILTER/PASS"};
    for (int r = 0; r < SIDELINES; r++) {
        fprintf(out_vcf, "##FILTER=<ID=%s,Description=\"%s\">\n",
                sideline_strs[r].data(), sideline_descs[r].data());
        declared.insert("FILTER/" + sideline_strs[r]);
    }

    // Declare every FILTER, INFO, and FORMAT field carried over from an input. A record's
    // site-level columns come from whichever callset owns it, so both callsets contribute; where
    // the two disagree about an ID, the query's declaration wins, as it does for a matched record.
    for (callset_t c : EnumRange<callset_t, CALLSET_SLOTS>{}) {
        if (this->callset_src_recs[c] == nullptr) continue;
        const srcRecords & src = *this->callset_src_recs[c];
        for (size_t i = 0; i < src.hdr_lines.size(); i++) {
            if (!declared.insert(src.hdr_keys[i]).second) continue;
            fprintf(out_vcf, "%s\n", src.hdr_lines[i].data());
        }
    }

    // The per-haplotype fields carry one value per allele of the sample's GT, which is what VCF
    // 4.4's Number=P declares. BCF_VL_P only reaches htslib in 1.23, so a consumer on any older
    // bcftools or pysam would report a cardinality error; Number=. produces byte-identical records
    // and merely gives up the declared cardinality, so the count and order are stated here instead.
    const char* per_allele = " One value per allele of this sample's GT, in GT allele order, "
            "'.' for a reference allele.";
    fprintf(out_vcf, "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"GenoType\">\n");
    fprintf(out_vcf, "##FORMAT=<ID=BD,Number=.,Type=String,Description=\"Benchmark Decision for call (TP/FP/FN, or N for a call that was not assessed).%s A record retained without evaluation was not assessed as a whole, so it carries a single N.\">\n", per_allele);
    fprintf(out_vcf, "##FORMAT=<ID=BC,Number=.,Type=Float,Description=\"Benchmark Credit (on the interval [0,1], based on sync group edit distance).%s\">\n", per_allele);
    fprintf(out_vcf, "##FORMAT=<ID=RD,Number=.,Type=Integer,Description=\"Reference edit Distance from truth within current sync group.%s\">\n", per_allele);
    fprintf(out_vcf, "##FORMAT=<ID=QD,Number=.,Type=Integer,Description=\"Query edit Distance from truth within current sync group.%s\">\n", per_allele);
    fprintf(out_vcf, "##FORMAT=<ID=BK,Number=.,Type=String,Description=\"BenchmarK category ('gm' if credit == 1, 'lm' if credit > 0, else '.').%s\">\n", per_allele);
    fprintf(out_vcf, "##FORMAT=<ID=QQ,Number=1,Type=Float,Description=\"variant Quality\">\n");
    fprintf(out_vcf, "##FORMAT=<ID=SC,Number=1,Type=Integer,Description=\"SuperCluster (index in contig)\">\n");
    fprintf(out_vcf, "##FORMAT=<ID=SG,Number=.,Type=Integer,Description=\"Sync Group (index in supercluster, for credit assignment).%s\">\n", per_allele);
    fprintf(out_vcf, "##FORMAT=<ID=PS,Number=1,Type=Integer,Description=\"Phase Set identifier (input, per-variant)\">\n");
    fprintf(out_vcf, "##FORMAT=<ID=PB,Number=1,Type=Integer,Description=\"Phase Block (output, per-supercluster, index in contig)\">\n");
    fprintf(out_vcf, "##FORMAT=<ID=BS,Number=1,Type=Integer,Description=\"Block Phase: 0 = PHASE_KEEP, 1 = PHASE_SWAP)\">\n");
    fprintf(out_vcf, "##FORMAT=<ID=VP,Number=1,Type=Integer,Description=\"Variant Phase: 0 = PHASE_ORIG, 1 = PHASE_SWAP, . = PHASE_NONE)\">\n");
    fprintf(out_vcf, "##FORMAT=<ID=FE,Number=1,Type=Integer,Description=\"Flip Error (a per-supercluster error)\">\n");
    fprintf(out_vcf, "##FORMAT=<ID=GE,Number=1,Type=String,Description=\"Genotype Error ('+' if 0/1 truth -> 1/1 query, '-' if 1/1 truth -> 0/1 query, '.' otherwise)\">\n");
    fprintf(out_vcf, "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tTRUTH\tQUERY\n");

    // write variants
    for (std::string ctg : this->contigs) {
        EnumArray<callset_t, int, CALLSET_SLOTS> ptrs{};
        EnumArray<callset_t, int, CALLSET_SLOTS> poss{};
        EnumArray<callset_t, bool, CALLSET_SLOTS> next{};
        std::shared_ptr<ctgPhaseblocks> ctg_pbs = this->phase_blocks[ctg];
        std::shared_ptr<ctgSuperclusters> ctg_scs = ctg_pbs->ctg_superclusters;
        auto & vars = ctg_pbs->ctg_superclusters->callset_vars;
        std::shared_ptr<ctgVariants> qvars = ctg_pbs->ctg_superclusters->callset_vars[QUERY];
        std::shared_ptr<ctgVariants> tvars = ctg_pbs->ctg_superclusters->callset_vars[TRUTH];

        // retained records are interleaved into the evaluated ones by position, not appended
        const ctg_sidelines side = sidelines_of(this->callset_sidelined, ctg);
        std::vector<int> side_ptrs = std::vector<int>(CALLSETS, 0);

        // flip/swap state comes from the query; these defaults hold on a contig it never calls on
        int phase_block = 0;
        phase_t block_state = PHASE_ORIG;
        bool flip_error = false;

        while ( ptrs[QUERY] < qvars->n || ptrs[TRUTH] < tvars->n) {

            // get next positions
            for (callset_t c : EnumRange<callset_t, CALLSET_SLOTS>{}) {
                poss[c] = ptrs[c] < int(vars[c]->poss.size()) ? 
                        vars[c]->poss[ptrs[c]] : 
                        std::numeric_limits<int>::max();
                if (ptrs[c] < int(vars[c]->types.size()) &&
                        (vars[c]->types[ptrs[c]] == TYPE_INS ||
                        vars[c]->types[ptrs[c]] == TYPE_DEL)) poss[c]--;
            }

            // set flags for next haps
            int pos = std::min(poss[QUERY], poss[TRUTH]);
            for (callset_t c : EnumRange<callset_t, CALLSET_SLOTS>{}) {
                next[c] = (poss[c] == pos);
            }
            write_sidelined(out_vcf, ctg, side, side_ptrs, pos);

            // update phasing
            if (ptrs[QUERY] < qvars->n) {
                phase_t phase = qvars->phases[ptrs[QUERY]];
                block_state = qvars->pb_phases[ptrs[QUERY]];

                // update switch/flip status
                if (block_state == PHASE_SWAP) {
                    if (phase == PHASE_ORIG) {
                        flip_error = true;
                    } else { // PHASE_SWAP or PHASE_NONE
                        flip_error = false;
                    }
                } else { // block_state == PHASE_ORIG
                    if (phase == PHASE_SWAP) {
                        flip_error = true;
                    } else { // PHASE_ORIG or PHASE_NONE
                        flip_error = false;
                    }
                }
            }

            // update supercluster and phase block
            int sc_idx = next[QUERY] ? vars[QUERY]->superclusters[ptrs[QUERY]] :
                                       vars[TRUTH]->superclusters[ptrs[TRUTH]];
            // the final entry is a past-the-end index, so stop advancing once it has been reached
            if (phase_block+1 < int(ctg_pbs->phase_blocks.size()) &&
                    ptrs[QUERY] >= ctg_pbs->phase_blocks[phase_block+1])
                phase_block++;

            /* if (next[QUERY] && ptrs[QUERY] < qvars->n) { */
            /*     fprintf(out_vcf, "orig_gt: %s\tmatched_gt: %s\tcredit: %.2f|%.2f\tref_dist: %d|%d\tphase: %s\n", */ 
            /*             gt_strs[vars[QUERY]->orig_gts[ptrs[QUERY]]].data(), */
            /*             gt_strs[vars[QUERY]->matched_gts[ptrs[QUERY]]].data(), */
            /*             vars[QUERY]->credit[HAP1][ptrs[QUERY]], */
            /*             vars[QUERY]->credit[HAP2][ptrs[QUERY]], */
            /*             vars[QUERY]->ref_ed[HAP1][ptrs[QUERY]], */
            /*             vars[QUERY]->ref_ed[HAP2][ptrs[QUERY]], */
            /*             phase_strs[vars[QUERY]->phases[ptrs[QUERY]]].data() */
            /*             ); */
            /* } */
            /* if (next[TRUTH] && ptrs[TRUTH] < tvars->n) { */
            /*     fprintf(out_vcf, "orig_gt: %s\n", */ 
            /*             gt_strs[vars[TRUTH]->orig_gts[ptrs[TRUTH]]].data()); */
            /* } */

            if (next[QUERY]) {
                // a positional tie between differing alleles is not a match: the two are written as
                // co-located records, the truth one on the next pass
                bool matched = next[TRUTH] &&
                        vars[QUERY]->refs[ptrs[QUERY]] == vars[TRUTH]->refs[ptrs[TRUTH]] &&
                        vars[QUERY]->alts[ptrs[QUERY]] == vars[TRUTH]->alts[ptrs[TRUTH]];

                // the query owns every record it appears on, so its source record supplies the
                // appended FORMAT keys and the truth sample has no values of its own for them
                const std::string pad = dot_fields(vars[QUERY]->src_fmt_keys(ptrs[QUERY]));
                vars[QUERY]->print_var_info(out_vcf, this->ref, ctg, ptrs[QUERY]);
                if (matched) {
                    vars[TRUTH]->print_var_sample(out_vcf, ptrs[TRUTH],
                            sc_idx, phase_block, block_state == PHASE_SWAP, flip_error, false, pad);
                } else {
                    vars[TRUTH]->print_var_empty(out_vcf, sc_idx, phase_block, false, pad);
                }
                vars[QUERY]->print_var_sample(out_vcf, ptrs[QUERY],
                        sc_idx, phase_block, block_state == PHASE_SWAP, flip_error, true,
                        vars[QUERY]->src_fmt_vals(ptrs[QUERY]));
                ptrs[QUERY]++;
                if (matched) ptrs[TRUTH]++;
            } else if (next[TRUTH]) {
                const std::string pad = dot_fields(vars[TRUTH]->src_fmt_keys(ptrs[TRUTH]));
                vars[TRUTH]->print_var_info(out_vcf, this->ref, ctg, ptrs[TRUTH]);
                vars[TRUTH]->print_var_sample(out_vcf, ptrs[TRUTH],
                        sc_idx, phase_block, block_state == PHASE_SWAP, flip_error, false,
                        vars[TRUTH]->src_fmt_vals(ptrs[TRUTH]));
                vars[QUERY]->print_var_empty(out_vcf, sc_idx, phase_block, true, pad);
                ptrs[TRUTH]++;
            } else {
                ERROR("No variants are selected next.");
            }
        }
        write_sidelined(out_vcf, ctg, side, side_ptrs, std::numeric_limits<int>::max());
    }

    // a contig absent from the BED file was dropped before evaluation, so the walk above never
    // reached it; nothing of it was evaluated, leaving its retained records the whole of its output
    for (const std::shared_ptr<ctgSideline> & only : retained_only) {
        std::vector<int> side_ptrs = std::vector<int>(CALLSETS, 0);
        write_sidelined(out_vcf, only->ctg, sidelines_of(this->callset_sidelined, only->ctg),
                side_ptrs, std::numeric_limits<int>::max());
    }
    fclose(out_vcf);
}


/**************************************************************************************************/


/**
 * @brief Constructs phaseblock container from supercluster data and runs phasing pipeline.
 * @param[in] clusterdata_ptr Supercluster data with contigs, lengths, and variants
 */
phaseblockData::phaseblockData(std::shared_ptr<superclusterData> clusterdata_ptr)
{
    // copy contigs and reference
    for (int i = 0; i < int(clusterdata_ptr->contigs.size()); i++) {
        std::string ctg = clusterdata_ptr->contigs[i];
        this->contigs.push_back(ctg);
        this->lengths.push_back(clusterdata_ptr->lengths[i]);
        this->phase_blocks[ctg] = std::shared_ptr<ctgPhaseblocks>(new ctgPhaseblocks());
    }
    this->ref = clusterdata_ptr->ref;
    this->callset_src_recs = clusterdata_ptr->callset_src_recs;
    this->callset_sidelined = clusterdata_ptr->callset_sidelined;

    // add pointers to superclusters
    for (const std::string & ctg : this->contigs) {
        this->phase_blocks[ctg]->ctg_superclusters = clusterdata_ptr->superclusters[ctg];
    }

    // fill in unset PS tags
    this->fix_phase_set_tags();

    // add phase blocks based on phase sets for each contig, after the unset tags are filled in: a
    // variant that declared no phase set of its own belongs to the block around it, and splitting
    // on its zero would report one block per such variant
    for (const std::string & ctg : this->contigs) {
        std::shared_ptr<ctgPhaseblocks> ctg_pbs = this->phase_blocks[ctg];
        std::shared_ptr<ctgVariants> qvars = ctg_pbs->ctg_superclusters->callset_vars[QUERY];
        int curr_phase_set = -1;
        for (int qvar_idx = 0; qvar_idx < qvars->n; qvar_idx++) {
            if (qvars->phase_sets[qvar_idx] != curr_phase_set) {
                ctg_pbs->phase_blocks.push_back(qvar_idx);
                ctg_pbs->n++;
                curr_phase_set = qvars->phase_sets[qvar_idx];
            }
        }
        ctg_pbs->phase_blocks.push_back(qvars->n);
    }

    // calculate phasings, flip, and switch errors
    this->phase();
    // calculate and fix allele count errors
    this->fix_allele_counts();
}


/**************************************************************************************************/


/**
 * @brief Propagates phase set tags to unphased and homozygous variants.
 * @note Must run before the phase block scan, phase(), and fix_allele_counts().
 * @todo Only set phase sets for 1|1 variants when unphased evaluation is added.
 */
void phaseblockData::fix_phase_set_tags() {

    for (callset_t ci : EnumRange<callset_t, CALLSET_SLOTS>{}) {
        // one span is recorded per phase set, so its count is the phase set count
        std::vector<int> phase_set_sizes;

        for (size_t ctg_idx = 0; ctg_idx < this->contigs.size(); ctg_idx++) { // for each contig
            const std::string & ctg = this->contigs[ctg_idx];
            std::shared_ptr<ctgSuperclusters> ctg_scs = this->phase_blocks[ctg]->ctg_superclusters;
            std::shared_ptr<ctgVariants> vars = ctg_scs->callset_vars[ci];

            // get first phase set (to backfill all preceding zeros)
            int first_phase_set = 0;
            for (int vi = 0; vi < vars->n; vi++) {
                if (vars->phase_sets[vi] != 0) {
                    first_phase_set = vars->phase_sets[vi];
                    break;
                }
            }

            // exit early if this contig has no phase sets
            if (first_phase_set == 0) {
                phase_set_sizes.push_back(this->lengths[ctg_idx]);
                continue;
            } else { // set phase set up until first PS
                for (int vi = 0; vi < vars->n; vi++) {
                    if (vars->phase_sets[vi] != 0) {
                        break;
                    }
                    vars->phase_sets[vi] = first_phase_set;
                }
            }

            int ps_beg = 0; int ps_end = 0;
            int phase_set = 0;
            for (int vi = 0; vi < vars->n; vi++) {
                if (vars->phase_sets[vi] != 0) { // variant is phased
                    if (vars->phase_sets[vi] != phase_set) { // new phase set, save old
                        if (phase_set) phase_set_sizes.push_back(ps_end - ps_beg);
                        phase_set = vars->phase_sets[vi];
                        ps_beg = vars->poss[vi];
                        ps_end = vars->poss[vi] + vars->rlens[vi];
                    } else { // same
                        ps_end = std::max(ps_end, vars->poss[vi] + vars->rlens[vi]);
                    }
                } else { // set unphased variant phase set to current phase set
                    // TODO: only set phase sets for 1|1 variants (for when we add unphased eval)
                    vars->phase_sets[vi] = phase_set;
                }
            }

            // add final phase set on contig
            phase_set_sizes.push_back(ps_end - ps_beg);
        }

        // calculate phaseset NG50
        size_t total_bases = 0;
        for (size_t i = 0; i < this->contigs.size(); i++) {
            total_bases += lengths[i];
        }

        int pb_ng50 = calc_ng50(phase_set_sizes, total_bases);

        if (g.verbosity >= 1) INFO("               %s phase sets: %zu",
                callset_strs[ci].data(), phase_set_sizes.size());
        if (g.verbosity >= 1) INFO("         %s phase block NG50: %d", 
                callset_strs[ci].data(), pb_ng50);
        if (g.verbosity >= 1) INFO("              %s total bases: %zu", 
                callset_strs[ci].data(), total_bases);
    }
}


/**************************************************************************************************/


/**
 * @brief Corrects matched genotypes to preserve allele counts matching original calls.
 * @note Records each variant's allele count error type on both callsets, and tracks and reports
 *       genotype error statistics (0/0->0/1, 1/1->0/1, etc.)
 * @throws ERROR if a query or truth variant's allele count error type is AC_UNKNOWN
 * @throws ERROR if the output genotype error TSV file cannot be opened for writing
 */
void phaseblockData::fix_allele_counts() {
    EnumArray<ac_errtype_t, EnumArray<sizeclass_t, int, SIZECLASS_SLOTS>,
            AC_ERRTYPE_SLOTS> allele_error_counts{};
    for (const std::string & ctg : this->contigs) {
        std::shared_ptr<ctgVariants> qvars = 
            this->phase_blocks[ctg]->ctg_superclusters->callset_vars[QUERY];

        for (int vi = 0; vi < qvars->n; vi++) {
            ac_errtype_t allele_count_errtype = qvars->set_allele_errtype(vi, true);
            if (allele_count_errtype == AC_UNKNOWN) {
                ERROR("Unknown variant allele count at %s:%d, %s -> %s", ctg.data(), qvars->poss[vi],
                        gt_strs[qvars->matched_gts[vi]].data(), gt_strs[qvars->orig_gts[vi]].data());
            }
            sizeclass_t vartype = qvars->get_vartype(vi);
            allele_error_counts[allele_count_errtype][vartype]++;
            allele_error_counts[allele_count_errtype][VARTYPE_ALL]++;

            // force 1|1 query variants to be evaluated as such
            if (qvars->orig_gts[vi] == GT_ALT_ALT) {
                qvars->matched_gts[vi] = qvars->orig_gts[vi];

                // if we're in a PHASE_SWAP phase block, we should swap data here since otherwise
                // there's no way based on matched_gt 1|1 to know to look at data from other hap
                if (qvars->pb_phases[vi] == PHASE_SWAP) {
                    std::swap(qvars->errtypes[HAP1][vi], qvars->errtypes[HAP2][vi]);
                    std::swap(qvars->sync_group[HAP1][vi], qvars->sync_group[HAP2][vi]);
                    std::swap(qvars->callq[HAP1][vi], qvars->callq[HAP2][vi]);
                    std::swap(qvars->ref_ed[HAP1][vi], qvars->ref_ed[HAP2][vi]);
                    std::swap(qvars->query_ed[HAP1][vi], qvars->query_ed[HAP2][vi]);
                    std::swap(qvars->credit[HAP1][vi], qvars->credit[HAP2][vi]);
                }
            }

            // force original 0|1 and 1|0 query variants to be evaluated as such, though truth differs
            // (matched_gt has allele count 2, orig_gt has allele count 1)
            else if (allele_count_errtype == AC_ERR_2_TO_1) {

                // for called 1|1 variants, keep variant call with better calculated credit
                if (qvars->credit[HAP1][vi] > qvars->credit[HAP2][vi]) {
                    qvars->set_var_matched_gt_on_hap(vi, HAP2, false);
                } else if (qvars->credit[HAP1][vi] < qvars->credit[HAP2][vi]) {
                    qvars->set_var_matched_gt_on_hap(vi, HAP1, false);
                } else { // default to current phasing
                    if (qvars->pb_phases[vi] == PHASE_ORIG) {
                        qvars->matched_gts[vi] = qvars->orig_gts[vi];
                    } else { // PHASE_SWAP
                        qvars->matched_gts[vi] = (qvars->orig_gts[vi] == GT_REF_ALT) ?
                            GT_ALT_REF : GT_REF_ALT;
                    }
                }
            // (matched_gt has allele count 0, orig_gt has allele count 1)
            } else if (allele_count_errtype == AC_ERR_0_TO_1) {
                // try to use hap with max credit
                if (qvars->credit[HAP1][vi] > qvars->credit[HAP2][vi]) {
                    qvars->set_var_matched_gt_on_hap(vi, HAP1, true);
                } else if (qvars->credit[HAP1][vi] < qvars->credit[HAP2][vi]) {
                    qvars->set_var_matched_gt_on_hap(vi, HAP2, true);
                } else { // default to current phasing
                    if (qvars->pb_phases[vi] == PHASE_ORIG) {
                        qvars->matched_gts[vi] = qvars->orig_gts[vi];
                    } else { // PHASE_SWAP
                        qvars->matched_gts[vi] = (qvars->orig_gts[vi] == GT_REF_ALT) ?
                            GT_ALT_REF : GT_REF_ALT;
                    }
                }
            }
        }

        // false negative errors can only be counted from the truth VCF
        std::shared_ptr<ctgVariants> tvars = 
            this->phase_blocks[ctg]->ctg_superclusters->callset_vars[TRUTH];
        for (int vi = 0; vi < tvars->n; vi++) {

            // the same value the query loop above records, read from the other side: a truth
            // record's own orig_gt supplies the truth allele count and its alignment-recovered
            // matched_gt the query's. This tallies nothing; the branches below own the summary.
            if (tvars->set_allele_errtype(vi, false) == AC_UNKNOWN) {
                ERROR("Unknown variant allele count at %s:%d, %s -> %s", ctg.data(), tvars->poss[vi],
                        gt_strs[tvars->orig_gts[vi]].data(), gt_strs[tvars->matched_gts[vi]].data());
            }

            sizeclass_t vartype = tvars->get_vartype(vi);
            if (tvars->orig_gts[vi] == GT_ALT_ALT) {
                if (tvars->errtypes[HAP1][vi] == ERRTYPE_FN && 
                        tvars->errtypes[HAP2][vi] == ERRTYPE_FN) {
                    allele_error_counts[AC_ERR_2_TO_0][vartype]++;
                    allele_error_counts[AC_ERR_2_TO_0][VARTYPE_ALL]++;
                }
            } else if (tvars->orig_gts[vi] == GT_REF_ALT) {
                if (tvars->errtypes[HAP2][vi] == ERRTYPE_FN) {
                    allele_error_counts[AC_ERR_1_TO_0][vartype]++;
                    allele_error_counts[AC_ERR_1_TO_0][VARTYPE_ALL]++;
                }
            } else if (tvars->orig_gts[vi] == GT_ALT_REF) {
                if (tvars->errtypes[HAP1][vi] == ERRTYPE_FN) {
                    allele_error_counts[AC_ERR_1_TO_0][vartype]++;
                    allele_error_counts[AC_ERR_1_TO_0][VARTYPE_ALL]++;
                }
            }
        }
    }

    if (g.verbosity >= 1) INFO(" ");
    if (g.verbosity >= 1) INFO("  Genotype Error Summary:");
    if (g.verbosity >= 1) INFO("  Truth -> Query");
    if (g.verbosity >= 1) INFO("    0/0 -> 0/1: %-8d  1 FP", 
            allele_error_counts[AC_ERR_0_TO_1][VARTYPE_ALL]);
    if (g.verbosity >= 1) INFO("    0/0 -> 1/1: %-8d  2 FP", 
            allele_error_counts[AC_ERR_0_TO_2][VARTYPE_ALL]);
    if (g.verbosity >= 1) INFO("    0/1 -> 0/0: %-8d  1 FN", 
            allele_error_counts[AC_ERR_1_TO_0][VARTYPE_ALL]);
    if (g.verbosity >= 1) INFO("    0/1 -> 0/1: %-8d  1 TP", 
            allele_error_counts[AC_ERR_1_TO_1][VARTYPE_ALL]);
    if (g.verbosity >= 1) INFO("    0/1 -> 1/1: %-8d  1 TP, 1 FP (False Homozygous)", 
            allele_error_counts[AC_ERR_1_TO_2][VARTYPE_ALL]);
    if (g.verbosity >= 1) INFO("    1/1 -> 0/0: %-8d  2 FN", 
            allele_error_counts[AC_ERR_2_TO_0][VARTYPE_ALL]);
    if (g.verbosity >= 1) INFO("    1/1 -> 0/1: %-8d  1 TP, 1 FN (False Heterozygous)", 
            allele_error_counts[AC_ERR_2_TO_1][VARTYPE_ALL]);
    if (g.verbosity >= 1) INFO("    1/1 -> 1/1: %-8d  2 TP", 
            allele_error_counts[AC_ERR_2_TO_2][VARTYPE_ALL]);
    if (g.verbosity >= 1) INFO(" ");
    this->write_genotype_error_summary(allele_error_counts);
}


/**************************************************************************************************/


/**
 * @brief Writes allele count error cross-tabulation table to TSV file.
 * @param[in] allele_error_counts 2D array indexed as [allele_count_errtype][vartype]
 * @throws ERROR if the output genotype error TSV file cannot be opened for writing
 */
void phaseblockData::write_genotype_error_summary(
        const EnumArray<ac_errtype_t, EnumArray<sizeclass_t, int, SIZECLASS_SLOTS>,
            AC_ERRTYPE_SLOTS> & allele_error_counts) {
    std::string out_genotype_errors_fn = g.out_prefix + "genotype-errors.tsv";
    FILE* out_genotype_errors = 0;
    if (g.verbosity >= 1) INFO("  Writing genotype error results to '%s'", out_genotype_errors_fn.data());
    out_genotype_errors = fopen(out_genotype_errors_fn.data(), "w");
    if (out_genotype_errors == NULL) {
        ERROR("Failed to open genotype error TSV file '%s'", out_genotype_errors_fn.data());
    }
    fprintf(out_genotype_errors, "VAR_TYPE\tALLELE_COUNT_0_TO_1\tALLELE_COUNT_0_TO_2\tALLELE_COUNT_1_TO_0\tALLELE_COUNT_1_TO_1\tALLELE_COUNT_1_TO_2\tALLELE_COUNT_2_TO_0\tALLELE_COUNT_2_TO_1\tALLELE_COUNT_2_TO_2\n");
    for (sizeclass_t vartype : EnumRange<sizeclass_t, SIZECLASS_SLOTS>{}) {
        fprintf(out_genotype_errors, "%s\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\n",
                vartype_strs[vartype].data(),
                allele_error_counts[AC_ERR_0_TO_1][vartype],
                allele_error_counts[AC_ERR_0_TO_2][vartype],
                allele_error_counts[AC_ERR_1_TO_0][vartype],
                allele_error_counts[AC_ERR_1_TO_1][vartype],
                allele_error_counts[AC_ERR_1_TO_2][vartype],
                allele_error_counts[AC_ERR_2_TO_0][vartype],
                allele_error_counts[AC_ERR_2_TO_1][vartype],
                allele_error_counts[AC_ERR_2_TO_2][vartype]);
    }
    fclose(out_genotype_errors);
}


/**************************************************************************************************/


/**
 * @brief Uses dynamic programming to find optimal phasing and detect switch/flip errors per contig.
 * @note Results stored in qvars->pb_phases and error lists in ctgPhaseblocks. Switches at phase set boundaries incur no cost.
 */
void phaseblockData::phase()
{
    // phase each contig separately
    for (const std::string & ctg : this->contigs) {
        std::shared_ptr<ctgPhaseblocks> ctg_pbs = this->phase_blocks[ctg];
        std::shared_ptr<ctgVariants> qvars = ctg_pbs->ctg_superclusters->callset_vars[QUERY];
        EnumArray<phase_t, std::vector<int>, PHASE_SLOTS> mat{};
        EnumArray<phase_t, std::vector<phaseptr_t>, PHASE_SLOTS> ptrs{};
        for (std::vector<int> & row : mat) row = std::vector<int>(qvars->n+1);
        for (std::vector<phaseptr_t> & row : ptrs) row = std::vector<phaseptr_t>(qvars->n+1);

        // calculate phasings for each variant
        for (int i = 0; i < qvars->n; i++) {
            if ((qvars->orig_gts[i] == GT_ALT_REF && qvars->matched_gts[i] == GT_ALT_REF) || // same
                    (qvars->orig_gts[i] == GT_REF_ALT && qvars->matched_gts[i] == GT_REF_ALT)) {
                qvars->phases[i] = PHASE_ORIG;
            } else if ((qvars->orig_gts[i] == GT_ALT_REF && qvars->matched_gts[i] == GT_REF_ALT) || // diff
                    (qvars->orig_gts[i] == GT_REF_ALT && qvars->matched_gts[i] == GT_ALT_REF)) {
                qvars->phases[i] = PHASE_SWAP;
            } else {
                qvars->phases[i] = PHASE_NONE;
            }
        }

        // forward pass
        for (int i = 0; i < qvars->n; i++) {

            // determine costs (penalized if this phasing deemed incorrect)
            EnumArray<phase_t, int, PHASE_SLOTS> costs{};
            switch (qvars->phases[i]) {
                case PHASE_ORIG:
                    costs[PHASE_ORIG] = 0; 
                    costs[PHASE_SWAP] = 1; 
                    break;
                case PHASE_SWAP:
                    costs[PHASE_ORIG] = 1; 
                    costs[PHASE_SWAP] = 0; 
                    break;
                case PHASE_NONE:
                    costs[PHASE_ORIG] = 0;
                    costs[PHASE_SWAP] = 0;
                    break;
            }

            // mat[i] ptrs[i] refer to points directly before cluster i, and if switch occurs it
            // is between cluster i and i+1

            // no cost for phase switches if on border of phase set
            int cost_swap = 1;
            if (i < qvars->n-1 && qvars->phase_sets[i] != qvars->phase_sets[i+1]) {
                cost_swap = 0;
            }

            for (phase_t phase : {PHASE_ORIG, PHASE_SWAP}) {
                phase_t other = other_phase(phase);
                if (mat[phase][i] + costs[phase] < mat[other][i] + costs[other] + cost_swap) {
                    mat[phase][i+1] = mat[phase][i] + costs[phase];
                    ptrs[phase][i+1] = PHASE_PTR_KEEP;
                }
                else {
                    mat[phase][i+1] = mat[other][i] + costs[other] + cost_swap;
                    ptrs[phase][i+1] = PHASE_PTR_SWAP;
                }
            }
        }

        // backwards pass
        if (qvars->n > 0) { // skip empty contigs

            // determine starting phase
            phase_t phase = PHASE_ORIG;
            if (mat[PHASE_SWAP][qvars->n] < mat[PHASE_ORIG][qvars->n])
                phase = PHASE_SWAP;

            int i = qvars->n;
            while (i > 0) {
                if (ptrs[phase][i] == PHASE_PTR_SWAP) {

                    // not a switch error if between phase sets
                    if (qvars->phase_sets[i] == qvars->phase_sets[i-1]) {
                        ctg_pbs->switches.push_back(i);
                        ctg_pbs->nswitches++;
                    }

                    phase = other_phase(phase);
                } else if (ptrs[phase][i] == PHASE_PTR_KEEP) { // within phase block
                    if (qvars->phases[i-1] != PHASE_NONE && qvars->phases[i-1] != phase) {
                        ctg_pbs->flips.push_back(i-1);
                        ctg_pbs->nflips++;
                    }
                }
                i--;
                qvars->pb_phases[i] = phase;
            }
            std::reverse(ctg_pbs->flips.begin(), ctg_pbs->flips.end());
            std::reverse(ctg_pbs->switches.begin(), ctg_pbs->switches.end());
        }
    }
    

    // print
    int switch_errors = 0;
    int flip_errors = 0;
    int phase_blocks = 0;
    int variants = 0;
    if (g.verbosity >= 1) INFO(" ");
    if (g.verbosity >= 1) INFO("  Contigs:");
    int id = 0;
    for (const std::string & ctg : this->contigs) {
        std::shared_ptr<ctgPhaseblocks> ctg_pbs = this->phase_blocks[ctg];
        std::shared_ptr<ctgVariants> qvars = ctg_pbs->ctg_superclusters->callset_vars[QUERY];

        // print errors per contig
        if (g.verbosity >= 1) {
            INFO("    [%2d] %s: %d switch errors, %d flip errors, %d phase blocks", id, ctg.data(), 
                    ctg_pbs->nswitches, ctg_pbs->nflips, ctg_pbs->n);
        }
        variants += qvars->n;
        switch_errors += ctg_pbs->nswitches;
        flip_errors += ctg_pbs->nflips;
        phase_blocks += ctg_pbs->n;
        id++;
    }
    int ng50 = this->calculate_ng50(false, false);
    int s_ngc50 = this->calculate_ng50(true, false);
    int sf_ngc50 = this->calculate_ng50(true, true);

    if (g.verbosity >= 1) INFO(" ");
    if (g.verbosity >= 1) INFO("             Total  phase blocks: %d", phase_blocks);
    if (g.verbosity >= 1) INFO("             Total switch errors: %d", switch_errors);
    if (g.verbosity >= 1) INFO("             Total   flip errors: %d", flip_errors);
    if (g.verbosity >= 1 && variants) 
        INFO("               Switch error rate: %.6f%%", 100*switch_errors/float(variants));
    if (g.verbosity >= 1 && variants) 
        INFO("                 Flip error rate: %.6f%%", 100*flip_errors/float(variants));
    if (g.verbosity >= 1) INFO("               Phase block  NG50: %d", ng50);
    if (g.verbosity >= 1) INFO("  (switch)     Phase block NGC50: %d", s_ngc50);
    if (g.verbosity >= 1) INFO("  (switchflip) Phase block NGC50: %d", sf_ngc50);
    this->write_phasing_summary(phase_blocks, switch_errors, flip_errors, variants, ng50, s_ngc50, sf_ngc50);
}


/**************************************************************************************************/


/**
 * @brief Writes detected switch and flip error locations and classifications to TSV file.
 * @note Error types: SWITCH_ERR, FLIP, SWITCH_AND_FLIP
 * @throws ERROR if the output switchflip TSV file cannot be opened for writing
 */
void phaseblockData::write_switchflips() {

    std::string out_sf_fn = g.out_prefix + "switchflips.tsv";
    FILE* out_sf = 0;
    if (g.verbosity >= 1) INFO("  Writing switchflip results to '%s'", out_sf_fn.data());
    out_sf = fopen(out_sf_fn.data(), "w");
    if (out_sf == NULL) {
        ERROR("Failed to open switchflip TSV file '%s'", out_sf_fn.data());
    }
    fprintf(out_sf, "CONTIG\tSTART\tSTOP\tSWITCH_TYPE\tVARIANT\tPHASE_BLOCK\n");

    // get sizes of each correct phase block (split on flips, not just switch)
    for (std::string ctg: this->contigs) {
        
        std::shared_ptr<ctgPhaseblocks> ctg_pbs = this->phase_blocks[ctg];
        std::shared_ptr<ctgVariants> qvars = ctg_pbs->ctg_superclusters->callset_vars[QUERY];

        int switch_idx = 0;
        int flip_idx = 0;
        int pb_idx = 1;
        if (qvars->n == 0) continue;

        // init start of this correct block
        int vi = 0;
        int next_vi = qvars->n; // default to last
        int beg = 0; int end = 0;
        switchtype_t type = SWITCHTYPE_NONE;

        while (true) {

            // get type and supercluster of next switch/flip/phaseset
            if (pb_idx < ctg_pbs->n && ctg_pbs->phase_blocks[pb_idx] <= next_vi) {
                type = SWITCHTYPE_SWITCH;
                next_vi = ctg_pbs->phase_blocks[pb_idx];
            }
            if (switch_idx < ctg_pbs->nswitches && ctg_pbs->switches[switch_idx] <= next_vi) {
                type = SWITCHTYPE_SWITCH_ERR;
                next_vi = ctg_pbs->switches[switch_idx];
            }
            // check flip last because it will cause a second switch
            if (flip_idx < ctg_pbs->nflips && ctg_pbs->flips[flip_idx] <= next_vi) {
                if (type == SWITCHTYPE_SWITCH)
                    type = SWITCHTYPE_SWITCH_AND_FLIP;
                else
                    type = SWITCHTYPE_FLIP;
                next_vi = ctg_pbs->flips[flip_idx];
            }
            if (type == SWITCHTYPE_NONE) { // all out-of-bounds
                break;
            }
            if (next_vi < vi) ERROR("Next variant (%d) is not after current variant (%d) in write_switchflips()", next_vi, vi);


            // get block(s)
            if (type == SWITCHTYPE_FLIP || type == SWITCHTYPE_SWITCH_AND_FLIP) {
                // switch could have occurred anywhere after last phased supercluster
                int left = next_vi-1;
                while (left > 0 && qvars->phases[left] == PHASE_NONE)
                    left--;
                if (left >= 0) {
                    beg = qvars->poss[left] + qvars->rlens[left];
                    end = qvars->poss[next_vi];
                    fprintf(out_sf, "%s\t%d\t%d\t%s\t%d\t%d\n", ctg.data(), beg, end, 
                            switch_strs[SWITCHTYPE_FLIP_BEG].data(), next_vi, pb_idx-1);
                }

                // switch could have occurred anywhere before next phased supercluster
                int right = next_vi+1;
                while (right < qvars->n-1 && qvars->phases[right] == PHASE_NONE)
                    right++;
                if (right < qvars->n) {
                    beg = qvars->poss[next_vi] + qvars->rlens[next_vi];
                    end = qvars->poss[right];
                    fprintf(out_sf, "%s\t%d\t%d\t%s\t%d\t%d\n", ctg.data(), beg, end, 
                            switch_strs[SWITCHTYPE_FLIP_END].data(), next_vi, pb_idx-1);
                }
                flip_idx++;
                if (type == SWITCHTYPE_SWITCH_AND_FLIP) pb_idx++;

            } else if (type == SWITCHTYPE_SWITCH) {
                // end of phase block, don't print anything since not an error
                pb_idx++;

            } else if (type == SWITCHTYPE_SWITCH_ERR) {
                // expand left/right from in between these variants
                int left = next_vi-1;
                while (left > 0 && qvars->phases[left] == PHASE_NONE)
                    left--;
                int right = next_vi;
                while (right < qvars->n-1 && qvars->phases[right] == PHASE_NONE)
                    right++;
                if (left >= 0 && right < qvars->n) {
                    beg = qvars->poss[left] + qvars->rlens[left];
                    end = qvars->poss[right];
                    fprintf(out_sf, "%s\t%d\t%d\t%s\t%d\t%d\n", ctg.data(), beg, end, 
                            switch_strs[SWITCHTYPE_SWITCH_ERR].data(), next_vi, pb_idx-1);
                }
                switch_idx++;
            }

            vi = next_vi;
            next_vi = qvars->n; // reset to end
            type = SWITCHTYPE_NONE;
        }
    }
    fclose(out_sf);
}


/**************************************************************************************************/


/**
 * @brief Writes phasing summary statistics to TSV file.
 * @param[in] phase_blocks Total count of phase blocks across all contigs
 * @param[in] switch_errors Total count of phase switch errors
 * @param[in] flip_errors Total count of phase flip errors
 * @param[in] variants Total number of phased query variants
 * @param[in] ng50 NG50 of phase blocks without any error breaks
 * @param[in] s_ngc50 NGC50 of phase blocks broken on switch errors
 * @param[in] sf_ngc50 NGC50 of phase blocks broken on switch and flip errors
 * @throws ERROR if the output phasing summary TSV file cannot be opened for writing
 */
void phaseblockData::write_phasing_summary(int phase_blocks, int switch_errors,
        int flip_errors, int variants, int ng50, int s_ngc50, int sf_ngc50) {
    std::string out_phasing_summary_fn = g.out_prefix + "phasing-summary.tsv";
    if (g.verbosity >= 1) INFO(" ");
    if (g.verbosity >= 1) INFO("  Writing phasing summary to '%s'", 
            out_phasing_summary_fn.data());
    FILE* out_phasing_summary = fopen(out_phasing_summary_fn.data(), "w");
    if (out_phasing_summary == NULL) {
        ERROR("Failed to open phasing summary TSV file '%s'", out_phasing_summary_fn.data());
    }
    fprintf(out_phasing_summary,
            "PHASE_BLOCKS\tSWITCH_ERRORS\tFLIP_ERRORS\tSWITCH_ERROR_RATE\tFLIP_ERROR_RATE\t"
            "NG_50\tSWITCH_NGC50\tSWITCHFLIP_NGC50\n");
    float switch_error_rate = variants ? 100*switch_errors/float(variants) : 0;
    float flip_error_rate = variants ? 100*flip_errors/float(variants) : 0;
    fprintf(out_phasing_summary, "%d\t%d\t%d\t%.6f%%\t%.6f%%\t%d\t%d\t%d", phase_blocks,
            switch_errors, flip_errors, switch_error_rate,
            flip_error_rate, ng50, s_ngc50, sf_ngc50);
    fclose(out_phasing_summary);
}


/**************************************************************************************************/


/**
 * @brief Returns the reference span of each correctly-phased block on one contig.
 *
 * Walks the contig's breakpoints in ascending variant order, closing a block at the variant before
 * each break and opening the next at the variant after it. Phase set boundaries always break a
 * block, since starting a new phase set is not an error; switch and flip errors break one only
 * when requested. A flip breaks twice, excising the flipped variant into a block of its own; a flip
 * on the last variant ends the walk, since no variant remains to open a block after it.
 *
 * @param[in] ctg_pbs The contig's phase set boundary, switch error, and flip index vectors.
 * @param[in] qvars The contig's query variants, supplying block bounds via poss and rlens.
 * @param[in] break_on_switch If true, split blocks at switch errors.
 * @param[in] break_on_flip If true, split blocks at flip errors.
 * @return Reference span of each block, in ascending position order; empty if the contig holds no
 *   query variants.
 * @throws ERROR if the breakpoint indices are not in ascending order.
 */
std::vector<int> correct_block_sizes(const std::shared_ptr<ctgPhaseblocks> & ctg_pbs,
        const std::shared_ptr<ctgVariants> & qvars, bool break_on_switch, bool break_on_flip) {

    std::vector<int> correct_blocks;
    if (qvars->n == 0) return correct_blocks;

    int switch_idx = 0;
    int flip_idx = 0;
    int pb_idx = 1;

    // init start of this correct block
    int vi = 0;
    int next_vi = qvars->n; // default to last
    int beg = qvars->poss[0];
    int end = 0;
    switchtype_t type = SWITCHTYPE_NONE;

    while (true) {

        // get type and supercluster of next switch/flip/phaseset
        if (pb_idx < ctg_pbs->n && ctg_pbs->phase_blocks[pb_idx] <= next_vi) {
            type = SWITCHTYPE_SWITCH;
            next_vi = ctg_pbs->phase_blocks[pb_idx];
        }
        if (break_on_switch && switch_idx < ctg_pbs->nswitches && ctg_pbs->switches[switch_idx] <= next_vi) {
            type = SWITCHTYPE_SWITCH_ERR;
            next_vi = ctg_pbs->switches[switch_idx];
        }
        // check flip last (takes preference due to <=) because it can cause two breaks (before/after)
        // NOTE: it is possible for one supercluster to have both a switch (new PS) and flip
        if (break_on_flip && flip_idx < ctg_pbs->nflips && ctg_pbs->flips[flip_idx] <= next_vi) {
            if (type == SWITCHTYPE_SWITCH)
                type = SWITCHTYPE_SWITCH_AND_FLIP;
            else
                type = SWITCHTYPE_FLIP;
            next_vi = ctg_pbs->flips[flip_idx];
        }
        if (type == SWITCHTYPE_NONE) { // all out-of-bounds
            break;
        }
        if (next_vi < vi) ERROR("Next variant (%d) is not after current variant (%d) in correct_block_sizes()", next_vi, vi);


        // get block(s)
        if (type == SWITCHTYPE_FLIP || type == SWITCHTYPE_SWITCH_AND_FLIP) {
            end = qvars->poss[next_vi-1] + qvars->rlens[next_vi-1];
            correct_blocks.push_back(end-beg);
            beg = qvars->poss[next_vi];

            end = qvars->poss[next_vi] + qvars->rlens[next_vi];
            correct_blocks.push_back(end-beg);
            // a flip on the last variant leaves no variants to open the following block with
            if (next_vi+1 == qvars->n) return correct_blocks;
            beg = qvars->poss[next_vi+1];
            flip_idx++;
            if (type == SWITCHTYPE_SWITCH_AND_FLIP) pb_idx++;
        } else if (type == SWITCHTYPE_SWITCH) {
            end = qvars->poss[next_vi-1] + qvars->rlens[next_vi-1];
            correct_blocks.push_back(end-beg);
            beg = qvars->poss[next_vi];
            pb_idx++;
        } else if (type == SWITCHTYPE_SWITCH_ERR) {
            end = qvars->poss[next_vi-1] + qvars->rlens[next_vi-1];
            correct_blocks.push_back(end-beg);
            beg = qvars->poss[next_vi];
            switch_idx++;
        }

        vi = next_vi;
        next_vi = qvars->n; // reset to end
        type = SWITCHTYPE_NONE;
    }
    end = qvars->poss[qvars->n-1] + qvars->rlens[qvars->n-1];
    correct_blocks.push_back(end-beg);
    return correct_blocks;
}


/**************************************************************************************************/


/**
 * @brief Calculates NGC50 of phase blocks, optionally broken at switch or flip errors.
 * @param[in] break_on_switch If true, split blocks at switch errors
 * @param[in] break_on_flip If true, split blocks at flip errors
 * @return NGC50 value (block length at 50% cumulative length), or 0 if no blocks exist
 */
int phaseblockData::calculate_ng50(bool break_on_switch, bool break_on_flip) {

    // get total bases in genome
    size_t total_bases = 0;
    for (size_t i = 0; i < this->contigs.size(); i++) {
        total_bases += lengths[i];
    }

    // get sizes of each correct phase block (split on flips, not just switch)
    std::vector<int> correct_blocks;
    for (const std::string & ctg: this->contigs) {
        std::shared_ptr<ctgPhaseblocks> ctg_pbs = this->phase_blocks[ctg];
        std::vector<int> ctg_blocks = correct_block_sizes(ctg_pbs,
                ctg_pbs->ctg_superclusters->callset_vars[QUERY], break_on_switch, break_on_flip);
        correct_blocks.insert(correct_blocks.end(), ctg_blocks.begin(), ctg_blocks.end());
    }
    return calc_ng50(correct_blocks, total_bases);
}
