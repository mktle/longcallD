#include "wavefront/wavefront_align.h"
#include "alignment/cigar.h"
#include "seq.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include "abpoa.h"
#include "edlib.h"
#include "call_var_main.h"
#include "deknot.h"
#include "bam_utils.h"
#include "utils.h"
#include "math_utils.h"
#include "align.h"
#include "kmer.h"
#include <pthread.h>
#include <inttypes.h>
#include <ctype.h>
#include <limits.h>

extern int LONGCALLD_VERBOSE;

/******************************************
 * DeKnot local-assembly region statistics *
 ******************************************/
// outcome of one assembly attempt (a deknot_sv_rescue_core call that reached the assembler)
enum {
    DK_OUT_COMMIT1 = 0, // committed a single consensus
    DK_OUT_COMMIT2,     // committed two haplotype walks
    DK_OUT_NO_WALKS,    // assembler returned no walks
    DK_OUT_MANY_WALKS,  // assembler returned more than two walks and none aligned spanning the noisy region
    DK_OUT_NO_GAIN,     // walks did not improve on the existing consensus count
    DK_OUT_POOR_ALN,    // a walk failed to align to the expanded reference
    DK_OUT_HIT_EDGE,    // an indel run crossed the small window edge (retried with the wide flank)
    DK_OUT_NO_SUPPORT,  // trimmed walks had no supporting noisy-region reads
    DK_OUT_PAIR_NO_SV,  // paired pass: no walk carries the SV the breakpoint pair predicts
    DK_OUT_PAIR_NO_JUNC,
    DK_OUT_N
};
static const char *dk_out_names[DK_OUT_N] = {"commit_1cons", "commit_2cons", "no_walks", "gt2_walks_none_span", "no_gain", "poor_aln", "hit_edge", "no_read_support", "paired_no_predicted_sv", "paired_no_junction_support"};
// which pass of deknot_sv_rescue invoked the assembler
enum { DK_PASS_SMALL = 0, DK_PASS_WIDE, DK_PASS_WIDE_RETRY, DK_PASS_PAIRED, DK_PASS_DEPTH, DK_PASS_N };
static const char *dk_pass_names[DK_PASS_N] = {"small_flank", "wide_flank", "wide_retry", "paired", "depth"};
// how a result with more than two walks was reduced to at most two
enum { DK_SEL_NONE = 0, DK_SEL_PAIR, DK_SEL_SINGLE, DK_SEL_N };
static const char *dk_sel_names[DK_SEL_N] = {"none", "pair", "single"};

typedef struct {
    pthread_mutex_t lock; // chunks are processed in parallel (kt_for)
    // noisy regions handed to deknot_sv_rescue, split by the number of consensus sequences they already had
    int64_t n_regions, n_regions_cons0, n_regions_cons1;
    int64_t n_regions_rescued; // regions where local assembly raised n_cons
    // attempts that stopped before running the assembler
    int64_t n_skip_no_breakpoint, n_skip_window, n_skip_depth;
    // attempts that ran the assembler
    int64_t n_asm, n_asm_pass[DK_PASS_N], n_asm_subsampled;
    int64_t n_walks_hist[4]; // 0, 1, 2, >=3 walks returned
    int64_t n_many_walks_hist[3]; // 3, 4, >=5 walks returned
    // no-walk attempts split by cause: dk_assemble_walks returned NULL (graph construction failed)
    // or a walk list with zero entries (phasing emitted nothing)
    int64_t n_no_walks_null, n_no_walks_empty;
    double no_walks_null_time, no_walks_empty_time;
    // two-walk results dropped to one consensus because one walk collected no supporting reads,
    // split by the outcome they ended with after the demotion
    int64_t n_demoted, n_demoted_outcome[DK_OUT_N];
    // results with more than two walks, split by what the walk selection kept (a pair differing by a large
    // indel over the noisy region, or the single best-aligning walk), with the outcome they ended with
    int64_t n_sel[DK_SEL_N], n_sel_outcome[DK_SEL_N][DK_OUT_N];
    // two-walk results collapsed to one consensus because the walks were near-identical over the noisy region
    int64_t n_collapsed, n_collapsed_outcome[DK_OUT_N];
    int64_t n_outcome[DK_OUT_N];
    // paired-breakpoint pass: regions with a partner, attempts that reached the assembler, commits carrying the
    // predicted SV, attempts rejected by the SV gate, attempts rejected before assembly (window/depth)
    int64_t n_pair_regions, n_pair_asm, n_pair_commit_sv, n_pair_gate_reject, n_pair_pre_reject;
    int64_t n_depth_runs, n_depth_asm, n_depth_commit_sv, n_depth_gate_reject;
    int64_t n_pair_junc_reject, n_depth_junc_reject;
    int64_t n_wfa_capped;
    double outcome_time_sum[DK_OUT_N];  // wall seconds inside dk_assemble_walks, per outcome
    int64_t outcome_reads_sum[DK_OUT_N]; // reads given to the assembler, per outcome
    int64_t reg_len_sum, reg_len_min, reg_len_max; // noisy region length (bp)
    int64_t win_len_sum, win_len_min, win_len_max; // expanded assembly window length (bp)
    int64_t n_reads_sum, n_reads_min, n_reads_max; // reads given to the assembler
    double asm_time_sum, asm_time_max;             // wall seconds inside dk_assemble_walks
} deknot_stats_t;

static deknot_stats_t dk_stats = { PTHREAD_MUTEX_INITIALIZER };

static void dk_stats_count(int64_t *counter) {
    pthread_mutex_lock(&dk_stats.lock);
    (*counter)++;
    pthread_mutex_unlock(&dk_stats.lock);
}

static void dk_stats_add_region(int n_cons_existing, int ret) {
    pthread_mutex_lock(&dk_stats.lock);
    dk_stats.n_regions++;
    if (n_cons_existing == 0) dk_stats.n_regions_cons0++; else dk_stats.n_regions_cons1++;
    if (ret > n_cons_existing) dk_stats.n_regions_rescued++;
    pthread_mutex_unlock(&dk_stats.lock);
}

#define DK_STATS_MINMAX(field, v) do { \
    dk_stats.field##_sum += (v); \
    if (dk_stats.n_asm == 1 || (v) < dk_stats.field##_min) dk_stats.field##_min = (v); \
    if (dk_stats.n_asm == 1 || (v) > dk_stats.field##_max) dk_stats.field##_max = (v); \
} while (0)

// record one assembly attempt and, in verbose mode, print its per-region line
// walks_null: dk_assemble_walks returned NULL; demoted: a two-walk result was dropped to one consensus
// sel: what the >2-walk selection kept; collapsed: a near-identical walk pair was reduced to one consensus
static double dk_aln_secs = 0.0, dk_read_secs = 0.0;
static void dk_stats_add_attempt(const bam_chunk_t *chunk, hts_pos_t noisy_reg_beg, hts_pos_t noisy_reg_end,
                                 hts_pos_t exp_beg, hts_pos_t exp_end, int n_reads, int n_noisy_reads, int subsampled,
                                 int n_walks, int walks_null, int demoted, int sel, int collapsed, double secs, int pass, int outcome) {
    int64_t reg_len = noisy_reg_end - noisy_reg_beg + 1, win_len = exp_end - exp_beg + 1;
    pthread_mutex_lock(&dk_stats.lock);
    dk_stats.n_asm++;
    dk_stats.n_asm_pass[pass]++;
    if (subsampled) dk_stats.n_asm_subsampled++;
    dk_stats.n_walks_hist[n_walks < 3 ? n_walks : 3]++;
    if (n_walks >= 3) dk_stats.n_many_walks_hist[n_walks < 5 ? n_walks - 3 : 2]++;
    if (n_walks == 0) {
        if (walks_null) { dk_stats.n_no_walks_null++; dk_stats.no_walks_null_time += secs; }
        else { dk_stats.n_no_walks_empty++; dk_stats.no_walks_empty_time += secs; }
    }
    if (demoted) { dk_stats.n_demoted++; dk_stats.n_demoted_outcome[outcome]++; }
    if (sel != DK_SEL_NONE) { dk_stats.n_sel[sel]++; dk_stats.n_sel_outcome[sel][outcome]++; }
    if (collapsed) { dk_stats.n_collapsed++; dk_stats.n_collapsed_outcome[outcome]++; }
    dk_stats.n_outcome[outcome]++;
    if (pass == DK_PASS_PAIRED) {
        dk_stats.n_pair_asm++;
        if (outcome == DK_OUT_COMMIT1 || outcome == DK_OUT_COMMIT2) dk_stats.n_pair_commit_sv++;
        else if (outcome == DK_OUT_PAIR_NO_SV) dk_stats.n_pair_gate_reject++;
        else if (outcome == DK_OUT_PAIR_NO_JUNC) dk_stats.n_pair_junc_reject++;
    } else if (pass == DK_PASS_DEPTH) {
        dk_stats.n_depth_asm++;
        if (outcome == DK_OUT_COMMIT1 || outcome == DK_OUT_COMMIT2) dk_stats.n_depth_commit_sv++;
        else if (outcome == DK_OUT_PAIR_NO_SV) dk_stats.n_depth_gate_reject++;
        else if (outcome == DK_OUT_PAIR_NO_JUNC) dk_stats.n_depth_junc_reject++;
    }
    dk_stats.outcome_time_sum[outcome] += secs;
    dk_stats.outcome_reads_sum[outcome] += n_reads;
    DK_STATS_MINMAX(reg_len, reg_len);
    DK_STATS_MINMAX(win_len, win_len);
    DK_STATS_MINMAX(n_reads, (int64_t)n_reads);
    dk_stats.asm_time_sum += secs;
    if (secs > dk_stats.asm_time_max) dk_stats.asm_time_max = secs;
    pthread_mutex_unlock(&dk_stats.lock);
    if (LONGCALLD_VERBOSE >= 1)
        fprintf(stderr, "DeKnotAsm %s:%" PRIi64 "-%" PRIi64 " %" PRIi64 " pass=%s win=%" PRIi64 "-%" PRIi64 " (%" PRIi64 ") reads=%d%s noisy_reads=%d walks=%d%s%s%s%s%s time=%.2fs aln=%.2fs reads_aln=%.2fs outcome=%s\n",
                chunk->tname, noisy_reg_beg, noisy_reg_end, reg_len, dk_pass_names[pass], exp_beg, exp_end, win_len,
                n_reads, subsampled ? "(subsampled)" : "", n_noisy_reads, n_walks, (n_walks == 0 && walks_null) ? "(null)" : "",
                demoted ? " demoted" : "", sel != DK_SEL_NONE ? " sel=" : "", sel != DK_SEL_NONE ? dk_sel_names[sel] : "",
                collapsed ? " collapsed" : "", secs, dk_aln_secs, dk_read_secs, dk_out_names[outcome]);
    dk_aln_secs = 0.0; dk_read_secs = 0.0;
}

#undef DK_STATS_MINMAX

void deknot_print_stats(void) {
    deknot_stats_t *s = &dk_stats;
    pthread_mutex_lock(&s->lock);
    _err_info("DeKnot local assembly summary:\n");
    fprintf(stderr, "  candidate noisy regions:  %" PRIi64 " (no consensus: %" PRIi64 ", one consensus: %" PRIi64 "), rescued: %" PRIi64 "\n",
            s->n_regions, s->n_regions_cons0, s->n_regions_cons1, s->n_regions_rescued);
    fprintf(stderr, "  skipped before assembly:  no breakpoint signature: %" PRIi64 ", window empty/too long: %" PRIi64 ", depth out of range: %" PRIi64 "\n",
            s->n_skip_no_breakpoint, s->n_skip_window, s->n_skip_depth);
    fprintf(stderr, "  assembly attempts:        %" PRIi64 " (small flank: %" PRIi64 ", wide flank: %" PRIi64 ", wide retry: %" PRIi64 ", paired breakpoints: %" PRIi64 ", depth-triggered: %" PRIi64 "), read set subsampled: %" PRIi64 "\n",
            s->n_asm, s->n_asm_pass[DK_PASS_SMALL], s->n_asm_pass[DK_PASS_WIDE], s->n_asm_pass[DK_PASS_WIDE_RETRY], s->n_asm_pass[DK_PASS_PAIRED], s->n_asm_pass[DK_PASS_DEPTH], s->n_asm_subsampled);
    fprintf(stderr, "  paired breakpoints:       regions paired: %" PRIi64 ", attempts: %" PRIi64 ", committed with predicted SV: %" PRIi64 ", rejected by SV gate: %" PRIi64 ", rejected by junction gate: %" PRIi64 ", rejected before assembly: %" PRIi64 "; depth-triggered: runs: %" PRIi64 ", attempts: %" PRIi64 ", committed: %" PRIi64 ", rejected by SV gate: %" PRIi64 ", rejected by junction gate: %" PRIi64 "\n",
            s->n_pair_regions, s->n_pair_asm, s->n_pair_commit_sv, s->n_pair_gate_reject, s->n_pair_junc_reject, s->n_pair_pre_reject,
            s->n_depth_runs, s->n_depth_asm, s->n_depth_commit_sv, s->n_depth_gate_reject, s->n_depth_junc_reject);
    fprintf(stderr, "  paired-pass walks stopped at the WFA score cap: %" PRIi64 "\n", s->n_wfa_capped);
    if (s->n_asm > 0) {
        double n = (double)s->n_asm;
        fprintf(stderr, "  noisy region length:      min %" PRIi64 ", mean %.0f, max %" PRIi64 " bp\n", s->reg_len_min, s->reg_len_sum / n, s->reg_len_max);
        fprintf(stderr, "  assembly window length:   min %" PRIi64 ", mean %.0f, max %" PRIi64 " bp\n", s->win_len_min, s->win_len_sum / n, s->win_len_max);
        fprintf(stderr, "  reads per assembly:       min %" PRIi64 ", mean %.1f, max %" PRIi64 "\n", s->n_reads_min, s->n_reads_sum / n, s->n_reads_max);
        fprintf(stderr, "  walks returned:           0: %" PRIi64 ", 1: %" PRIi64 ", 2: %" PRIi64 ", >=3: %" PRIi64 "\n",
                s->n_walks_hist[0], s->n_walks_hist[1], s->n_walks_hist[2], s->n_walks_hist[3]);
        fprintf(stderr, "  outcomes:                ");
        for (int i = 0; i < DK_OUT_N; ++i) fprintf(stderr, " %s: %" PRIi64 "%s", dk_out_names[i], s->n_outcome[i], i + 1 < DK_OUT_N ? "," : "\n");
        fprintf(stderr, "  assembly time:            total %.1f s, mean %.2f s, max %.2f s\n", s->asm_time_sum, s->asm_time_sum / n, s->asm_time_max);
        // attribution of the failure buckets
        fprintf(stderr, "  no walks by cause:        assembler returned NULL: %" PRIi64 " (mean %.2f s), empty walk list: %" PRIi64 " (mean %.2f s)\n",
                s->n_no_walks_null, s->n_no_walks_null > 0 ? s->no_walks_null_time / (double)s->n_no_walks_null : 0.0,
                s->n_no_walks_empty, s->n_no_walks_empty > 0 ? s->no_walks_empty_time / (double)s->n_no_walks_empty : 0.0);
        fprintf(stderr, "  >=3 walks by count:       3: %" PRIi64 ", 4: %" PRIi64 ", >=5: %" PRIi64 "\n",
                s->n_many_walks_hist[0], s->n_many_walks_hist[1], s->n_many_walks_hist[2]);
        fprintf(stderr, "  two walks demoted to one: %" PRIi64 " (one walk had no read support)", s->n_demoted);
        for (int i = 0, first = 1; i < DK_OUT_N; ++i) {
            if (s->n_demoted_outcome[i] == 0) continue;
            fprintf(stderr, "%s %s: %" PRIi64, first ? ", ended as:" : ",", dk_out_names[i], s->n_demoted_outcome[i]);
            first = 0;
        }
        fprintf(stderr, "\n");
        fprintf(stderr, "  >=3 walks selection:      pair: %" PRIi64 ", single: %" PRIi64 ", none spanning the region: %" PRIi64 "\n",
                s->n_sel[DK_SEL_PAIR], s->n_sel[DK_SEL_SINGLE], s->n_outcome[DK_OUT_MANY_WALKS]);
        for (int k = DK_SEL_PAIR; k < DK_SEL_N; ++k) {
            if (s->n_sel[k] == 0) continue;
            fprintf(stderr, "    %s selected ended as:", dk_sel_names[k]);
            for (int i = 0, first = 1; i < DK_OUT_N; ++i) {
                if (s->n_sel_outcome[k][i] == 0) continue;
                fprintf(stderr, "%s %s: %" PRIi64, first ? "" : ",", dk_out_names[i], s->n_sel_outcome[k][i]);
                first = 0;
            }
            fprintf(stderr, "\n");
        }
        fprintf(stderr, "  identical pair collapsed: %" PRIi64 " (two walks near-identical over the noisy region, kept one)", s->n_collapsed);
        for (int i = 0, first = 1; i < DK_OUT_N; ++i) {
            if (s->n_collapsed_outcome[i] == 0) continue;
            fprintf(stderr, "%s %s: %" PRIi64, first ? ", ended as:" : ",", dk_out_names[i], s->n_collapsed_outcome[i]);
            first = 0;
        }
        fprintf(stderr, "\n");
        fprintf(stderr, "  cost per outcome:        ");
        for (int i = 0, first = 1; i < DK_OUT_N; ++i) {
            if (s->n_outcome[i] == 0) continue;
            double m = (double)s->n_outcome[i];
            fprintf(stderr, "%s %s: %" PRIi64 " (%.1f s total, mean %.2f s, mean %.1f reads)", first ? "" : ",", dk_out_names[i], s->n_outcome[i],
                    s->outcome_time_sum[i], s->outcome_time_sum[i] / m, s->outcome_reads_sum[i] / m);
            first = 0;
        }
        fprintf(stderr, "\n");
    }
    pthread_mutex_unlock(&s->lock);
}


// INS/DEL : left-most
// short   :    XXXXX                      [TSD]XXXXX
// long    :    XXXXX[TSD]----------[polyA][TSD]XXXXX
// or:
// short   :    XXXXX                      [TSD]XXXXX
// long    :    XXXXX[TSD][polyT]----------[TSD]XXXXX
// INS/DEL : right-most (not considered in this function)
// short   :    XXXXX[TSD]                      XXXXX
// long    :    XXXXX[TSD]----------[polyA][TSD]XXXXX
// or:
// short   :    XXXXX[TSD]                      XXXXX
// long    :    XXXXX[TSD][polyT]----------[TSD]XXXXX
int collect_te_info(const call_var_opt_t *opt, int var_type, uint8_t *gap_seq, uint8_t *flank_ref_seq, int gap_len, hts_pos_t gap_pos,
                    uint8_t **tsd_seq, hts_pos_t *tsd_pos1, hts_pos_t *tsd_pos2, int *tsd_polya_len, int *te_seq_i, int *te_is_rev) {
    *tsd_seq = NULL; *tsd_pos1 = -1; *tsd_pos2 = -1; *tsd_polya_len = -1; *te_seq_i = -1; *te_is_rev = 0;
    int min_tsd_len = opt->min_tsd_len, max_tsd_len = opt->max_tsd_len;
    int min_polya_len = opt->min_polya_len; float min_polya_ratio = opt->min_polya_ratio;
    // compare long_seq[gap_start_i:] with ref_seq[gap_end_i+1:]
    int tsd_len = 0, n_mis = 0, max_allow_mis = 1;
    for (int i = 0; i < gap_len; i++) {
        uint8_t base1 = gap_seq[i], base2 = flank_ref_seq[i];
        if (base1 == base2) tsd_len = i+1;
        else {
            n_mis++;
            if (n_mis > max_allow_mis) break;
        }
        if (tsd_len > max_tsd_len) break;
    }
    int has_tsd = 0, has_polya = 0, max_search_polya_len = 20;
    if (tsd_len >= min_tsd_len && tsd_len <= max_tsd_len) {
        has_tsd = 1;
        // look for polyA
        for (int polya_len = 0, polya = 0, i = gap_len-1; i >= 0; i--) {
            polya_len++;
            if (gap_seq[i] == 0) { polya++;
                if (polya_len >= min_polya_len && polya >= min_polya_ratio*polya_len) {
                    has_polya = 1;
                    *tsd_polya_len = polya_len; // update tsd_polya_len
                }
            } else if (polya_len > max_search_polya_len) break; // stop searching if too long
        }
        if (has_polya == 0) { // look for polyT
            for (int polyt_len = 0, polyt = 0, i = tsd_len; i < gap_len; i++) {
                polyt_len++;
                if (gap_seq[i] == 3) {
                    polyt++;
                    if (polyt_len >= min_polya_len && polyt >= min_polya_ratio*polyt_len) {
                        has_polya = 1;
                        *tsd_polya_len = -polyt_len; // update tsd_polya_len
                    }
                } else if (polyt_len > max_search_polya_len) break;
            }
        }
    }
    // check if INS/DEL is TE
    if (has_tsd && has_polya) {
        if (opt->n_te_seqs > 0) *te_seq_i = check_te_seq(opt, gap_seq, gap_len, te_is_rev);
        *tsd_seq = (uint8_t*)malloc(tsd_len*sizeof(uint8_t));
        for (int i = 0; i < tsd_len; i++) (*tsd_seq)[i] = flank_ref_seq[i];
        *tsd_pos1 = gap_pos; 
        if (var_type == BAM_CDEL) *tsd_pos2 = gap_pos+gap_len; else *tsd_pos2 = -1;
        return tsd_len;
    } else return 0;
}


// for somatic TE variant:
int collect_te_info_from_var(const call_var_opt_t *opt, bam_chunk_t *chunk, cand_var_t *var) {
    if (var->checked_tsd) return var->tsd_len; // already checked
    var->checked_tsd = 1; // mark as checked
    if (!((var->var_type == BAM_CINS && var->alt_len >= opt->min_sv_len) ||
        (var->var_type == BAM_CDEL && var->ref_len >= opt->min_sv_len)))
        return 0;
    // check if have TSD & polyA/T
    char *ref_seq = chunk->ref_seq; hts_pos_t ref_beg = chunk->ref_beg, ref_end = chunk->ref_end;
    uint8_t *gap_seq, *flank_ref_seq; int gap_len;
    if (var->var_type == BAM_CINS) {
        gap_len = var->alt_len;
        gap_seq = (uint8_t*)malloc(var->alt_len*sizeof(uint8_t));
        flank_ref_seq = (uint8_t*)malloc(var->alt_len*sizeof(uint8_t));
        for (int i = 0; i < gap_len; ++i) {
            gap_seq[i] = var->alt_seq[i];
            flank_ref_seq[i] = get_bseq1(ref_seq, ref_beg, ref_end, var->pos+i);
        }
    } else {
        gap_len = var->ref_len;
        gap_seq = (uint8_t*)malloc(var->ref_len*sizeof(uint8_t));
        flank_ref_seq = (uint8_t*)malloc(var->ref_len*sizeof(uint8_t));
        for (int i = 0; i < gap_len; ++i) {
            gap_seq[i] = get_bseq1(ref_seq, ref_beg, ref_end, var->pos+i);
            flank_ref_seq[i] = get_bseq1(ref_seq, ref_beg, ref_end, var->pos+i+gap_len);
        }
    }
    int tsd_len = 0, tsd_polya_len=0, te_seq_i=-1, te_is_rev=0;
    uint8_t *tsd_seq=NULL; hts_pos_t tsd_pos1=-1, tsd_pos2=-1; 
    tsd_len = collect_te_info(opt, var->var_type, gap_seq, flank_ref_seq, gap_len, var->pos,
                              &tsd_seq, &tsd_pos1, &tsd_pos2, &tsd_polya_len, &te_seq_i, &te_is_rev);
    free(gap_seq); free(flank_ref_seq);
    if (tsd_len > 0) { // copy tsd info to var
        var->tsd_seq = tsd_seq;
        var->tsd_len = tsd_len;
        var->polya_len = tsd_polya_len;
        var->tsd_pos1 = tsd_pos1;
        var->tsd_pos2 = tsd_pos2;
        var->te_seq_i = te_seq_i;
        var->te_is_rev = te_is_rev;
        return tsd_len;
    } else {
        var->tsd_seq = NULL; var->tsd_len = 0; var->tsd_pos1 = -1; var->tsd_pos2 = -1; var->polya_len = 0;
        var->te_seq_i = -1; var->te_is_rev = 0;
        return 0;
    }
}

// for germline TE variant:
// input  : cons' aln-str, or raw read' aln-str (mosaic variant)
// return : 0: not have TSD & polyA, 1: have TSD & polyA
//        tsd_len: length of TSD
//        tsd_target_pos1/2: position of 1st/2nd TSD on target
int collect_te_info_from_cons(const call_var_opt_t *opt, bam_chunk_t *chunk, hts_pos_t gap_ref_start, int msa_gap_start, int var_type, int gap_len, uint8_t *cons_msa_seq, 
                              uint8_t **tsd_seq, hts_pos_t *tsd_pos1, hts_pos_t *tsd_pos2, int *tsd_polya_len, int *te_seq_i, int *te_is_rev) {
    xassert(var_type == BAM_CINS || var_type == BAM_CDEL, "Var_type should be INS/DEL\n");
    // in case flanking sequence in ref_msa_seq is not long enough
    char *ref_seq = chunk->ref_seq; hts_pos_t ref_beg = chunk->ref_beg, ref_end = chunk->ref_end;
    uint8_t *gap_seq = (uint8_t*)malloc(gap_len*sizeof(uint8_t));
    uint8_t *flank_ref_seq = (uint8_t*)malloc(gap_len*sizeof(uint8_t));
    if (var_type == BAM_CINS) {
        for (int i = 0; i < gap_len; ++i) {
            gap_seq[i] = cons_msa_seq[msa_gap_start+i];
            flank_ref_seq[i] = get_bseq1(ref_seq, ref_beg, ref_end, gap_ref_start+i);
        }
    } else {
        for (int i = 0; i < gap_len; ++i) {
            gap_seq[i] = get_bseq1(ref_seq, ref_beg, ref_end, gap_ref_start+i);
            flank_ref_seq[i] = get_bseq1(ref_seq, ref_beg, ref_end, gap_ref_start+i+gap_len);
        }
    }
    int tsd_len = 0;
    tsd_len = collect_te_info(opt, var_type, gap_seq, flank_ref_seq, gap_len, gap_ref_start,
                    tsd_seq, tsd_pos1, tsd_pos2, tsd_polya_len, te_seq_i, te_is_rev);
    free(gap_seq); free(flank_ref_seq);
    return tsd_len;
}

int* edlibAlignmentToXID(const unsigned char* const alignment, const int alignmentLength) {
    // Maps move code from alignment to char in cigar.
    //                      0    1    2    3
    // moveCodeToChar[] : {'=', 'I', 'D', 'X'};
    int *xid = (int*) calloc(4, sizeof(int));

    int i;
    for (i = 0; i < alignmentLength; i++) {
        if (alignment[i] == EDLIB_EDOP_MATCH) {
            xid[0]++; // equal
        } else if (alignment[i] == EDLIB_EDOP_MISMATCH) {
            xid[1]++; // mismatch
        } else if (alignment[i] == EDLIB_EDOP_INSERT) {
            xid[2]++; // insertion
        } else if (alignment[i] == EDLIB_EDOP_DELETE) {
            xid[3]++; // deletion
        } else {
            free(xid);
            return NULL; // invalid alignment code
        }
    }
    return xid;
}

// count number of mismatches and gaps from edlib alignment
int edlibAlignmentToXGAPS(const unsigned char* const alignment, const int alignmentLength) {
    if (LONGCALLD_VERBOSE >= 3) {
        fprintf(stderr, "Edlib alignment: ");
        for (int i = 0; i < alignmentLength; i++) {
            fprintf(stderr, "%d", alignment[i]);
        }
        fprintf(stderr, "\n");
    }
    int n_gaps = 0, n_mismatch = 0;
    for (int i = 0; i < alignmentLength; i++) {
            if (alignment[i] == EDLIB_EDOP_MATCH) {
            continue;
        } else if (alignment[i] == EDLIB_EDOP_MISMATCH) {
            n_mismatch++;
        } else if (alignment[i] == EDLIB_EDOP_INSERT || alignment[i] == EDLIB_EDOP_DELETE) {
            if (i == 0 || alignment[i - 1] != alignment[i]) n_gaps++;
        }
    }
    return n_mismatch + n_gaps;
}

int edlib_edit_distance(uint8_t *target, int tlen, uint8_t *query, int qlen) {
    EdlibAlignResult result = edlibAlign((const char*)query, qlen, (const char*)target, tlen,
                                         edlibNewAlignConfig(-1, EDLIB_MODE_NW, EDLIB_TASK_DISTANCE, NULL, 0));
    if (result.status != EDLIB_STATUS_OK) {
        edlibFreeAlignResult(result);
        return -1;
    }
    int score = result.editDistance;
    edlibFreeAlignResult(result);
    return score;
}

int edlib_xgaps(uint8_t *target, int tlen, uint8_t *query, int qlen) {
    EdlibAlignResult result = edlibAlign((const char*)query, qlen, (const char*)target, tlen,
                                         edlibNewAlignConfig(-1, EDLIB_MODE_NW, EDLIB_TASK_PATH, NULL, 0));
    if (result.status != EDLIB_STATUS_OK) {
        edlibFreeAlignResult(result);
        return -1;
    }
    int n_x_gaps = edlibAlignmentToXGAPS(result.alignment, result.alignmentLength);
    edlibFreeAlignResult(result);
    return n_x_gaps;
}

int edlib_end2end_aln(uint8_t *target, int tlen, uint8_t *query, int qlen, int *n_eq, int *n_xid) {
    EdlibAlignResult result = edlibAlign((const char*)query, qlen, (const char*)target, tlen,
                                         edlibNewAlignConfig(-1, EDLIB_MODE_NW, EDLIB_TASK_PATH, NULL, 0));
    if (result.status != EDLIB_STATUS_OK) {
        edlibFreeAlignResult(result);
        return -1;
    }
    if (n_eq != NULL && n_xid != NULL) { // collect equal and XID bases
        int *xid = edlibAlignmentToXID(result.alignment, result.alignmentLength);
        if (xid != 0) {
            *n_eq = xid[0]; *n_xid = xid[1] + xid[2] + xid[3];
            free(xid);
        } else {
            *n_eq = -1; *n_xid = -1;
        }
    }
    int score = result.editDistance;
    edlibFreeAlignResult(result);
    return score;
}

// infix/HW alignment using edlib
int edlib_infix_aln(uint8_t *target, int tlen, uint8_t *query, int qlen, int *n_eq, int *n_xid) {
    EdlibAlignResult result = edlibAlign((const char*)query, qlen, (const char*)target, tlen,
                                         edlibNewAlignConfig(-1, EDLIB_MODE_HW, EDLIB_TASK_PATH, NULL, 0));
    if (result.status != EDLIB_STATUS_OK) {
        edlibFreeAlignResult(result);
        return -1;
    }
    if (n_eq != NULL && n_xid != NULL) { // collect equal and XID bases
        int *xid = edlibAlignmentToXID(result.alignment, result.alignmentLength);
        if (xid != 0) {
            *n_eq = xid[0]; *n_xid = xid[1] + xid[2] + xid[3];
            free(xid);
        } else {
            *n_eq = -1; *n_xid = -1;
        }
    }
    int score = result.editDistance;
    edlibFreeAlignResult(result);
    return score;
}

int wfa_collect_pretty_alignment(cigar_t* const cigar,
                                 const uint8_t* const pattern,
                                 const int pattern_length,
                                 const uint8_t* const text,
                                 const int text_length,
                                 uint8_t **_pattern_alg, uint8_t **_text_alg) {
    // Parameters
    char* const operations = cigar->operations;
    const int begin_offset = cigar->begin_offset;
    const int end_offset = cigar->end_offset;
    // Allocate alignment buffers
    const int max_buffer_length = text_length + pattern_length + 1;
    uint8_t *mem = (uint8_t*)calloc(2*max_buffer_length,1);
    uint8_t *pattern_alg = mem;
    uint8_t *text_alg = mem + max_buffer_length;
    // Compute alignment buffers
    int i, alg_pos = 0, pattern_pos = 0, text_pos = 0;
    for (i=begin_offset;i<end_offset;++i) {
        switch (operations[i]) {
            case 'M':
                if (pattern[pattern_pos] != text[text_pos]) {
                    pattern_alg[alg_pos] = pattern[pattern_pos];
                    text_alg[alg_pos++] = text[text_pos];
                } else {
                    pattern_alg[alg_pos] = pattern[pattern_pos];
                    text_alg[alg_pos++] = text[text_pos];
                }
                pattern_pos++; text_pos++;
                break;
            case 'X':
                if (pattern[pattern_pos] != text[text_pos]) {
                    pattern_alg[alg_pos] = pattern[pattern_pos++];
                    text_alg[alg_pos++] = text[text_pos++];
                } else {
                    pattern_alg[alg_pos] = pattern[pattern_pos++];
                    text_alg[alg_pos++] = text[text_pos++];
                }
                break;
            case 'I':
                pattern_alg[alg_pos] = 5; //'-';
                text_alg[alg_pos++] = text[text_pos++];
                break;
            case 'D':
                pattern_alg[alg_pos] = pattern[pattern_pos++];
                text_alg[alg_pos++] = 5; //'-';
                break;
            default: break;
        }
    }
    // assert(pattern_pos == pattern_length && text_pos == text_length);
    *_pattern_alg = pattern_alg; *_text_alg = text_alg;
    return alg_pos;
}

// return alignment score
int wfa_heuristic_aln(uint8_t *pattern, int plen, uint8_t *text, int tlen, 
                      int a, int b, int q, int e, int q2, int e2, int *n_eq, int *n_xid) {
    wavefront_aligner_attr_t attributes = wavefront_aligner_attr_default;
    attributes.memory_mode = wavefront_memory_ultralow;
    attributes.distance_metric = gap_affine_2p;
    attributes.affine2p_penalties.match = 0;
    attributes.affine2p_penalties.mismatch = b;
    attributes.affine2p_penalties.gap_opening1 = q;
    attributes.affine2p_penalties.gap_extension1 = e;
    attributes.affine2p_penalties.gap_opening2 = q2;
    attributes.affine2p_penalties.gap_extension2 = e2;
    attributes.alignment_scope = compute_alignment;
    attributes.alignment_form.span = alignment_end2end;
    // xdrops
    attributes.heuristic.strategy = wf_heuristic_xdrop;
    attributes.heuristic.xdrop = MAX_OF_TWO(100, (int)((plen + tlen) * 0.1)); 
    attributes.heuristic.steps_between_cutoffs = 100;
    // Initialize Wavefront Aligner
    wavefront_aligner_t* const wf_aligner = wavefront_aligner_new(&attributes);
    // Align
    wavefront_align(wf_aligner, (const char*)pattern, plen, (const char*)text, tlen);
    // collect score
    int score = wf_aligner->cigar->score; //cigar_score_gap_affine2p(wf_aligner->cigar, &attributes.affine2p_penalties);
    if (n_eq != NULL && n_xid != NULL) { // collect equal and XID bases
        uint32_t *cigar; int cigar_length;
        // collect cigar
        cigar_get_CIGAR(wf_aligner->cigar, true, &cigar, &cigar_length);
        *n_eq = 0, *n_xid = 0;
        for (int i = 0; i < cigar_length; ++i) {
            int op = cigar[i] & 0xf; int len = cigar[i] >> 4;
            if (op == BAM_CEQUAL) *n_eq += len;
            else if (op == BAM_CDIFF || op == BAM_CINS || op == BAM_CDEL) *n_xid += len;
        }
    }
    // Free
    wavefront_aligner_delete(wf_aligner); 
    return score;
}

// cons vs ref: ultra-low mem + dual-gap + no heuristic to ensure accurate alignment
// full read/cons vs full read: ultra-low mem + affine-gap + wf-adaptive to ensure accurate alignment and speed up
// full read/cons vs partial read: high mem (default) + affine-gap + heuristic (xdrop/zdrop) to speed up
static __thread int wfa_max_steps = INT_MAX;
static __thread int wfa_last_capped = 0;
int wfa_end2end_aln_mem(int mem_ultralow, int pattern_ends_free, uint8_t *pattern, int plen, uint8_t *text, int tlen,
                    int gap_aln, int b, int q, int e, int q2, int e2, int heuristic, int affine_gap, // heuristic: 0: no, 1: default, 2: zdrop
                    uint32_t **cigar_buf, int *cigar_length, uint8_t **pattern_alg, uint8_t **text_alg, int *alg_length) {
    // double realtime0 = realtime();
    // fprintf(stderr, "WFA-end2end %d vs %d\n", plen, tlen);
    wavefront_aligner_attr_t attributes = wavefront_aligner_attr_default;
    if (mem_ultralow) attributes.memory_mode = wavefront_memory_ultralow;
    if (affine_gap == LONGCALLD_WFA_AFFINE_2P) {
        attributes.distance_metric = gap_affine_2p;
        attributes.affine2p_penalties.match = 0; // -a;
        attributes.affine2p_penalties.mismatch = b;
        attributes.affine2p_penalties.gap_opening1 = q;
        attributes.affine2p_penalties.gap_extension1 = e;
        attributes.affine2p_penalties.gap_opening2 = q2;
        attributes.affine2p_penalties.gap_extension2 = e2;
    } else {
        attributes.distance_metric = gap_affine;
        attributes.affine_penalties.match = 0; // -a;
        attributes.affine_penalties.mismatch = b;
        attributes.affine_penalties.gap_opening = q;
        attributes.affine_penalties.gap_extension = e;
    }
    attributes.alignment_scope = compute_alignment;
    if (pattern_ends_free > 0) {
        attributes.alignment_form.span = alignment_endsfree;
        attributes.alignment_form.pattern_begin_free = pattern_ends_free;
        attributes.alignment_form.pattern_end_free = pattern_ends_free;
        attributes.alignment_form.text_begin_free = 0;
        attributes.alignment_form.text_end_free = 0;
    } else attributes.alignment_form.span = alignment_end2end;
    if (heuristic == LONGCALLD_WFA_NO_HEURISTIC) 
        attributes.heuristic.strategy = wf_heuristic_none;
    if (heuristic == LONGCALLD_WFA_ADAPTIVE) // default heuristic
        attributes.heuristic.strategy = wf_heuristic_wfadaptive;
    else if (heuristic == LONGCALLD_WFA_ZDROP) { // Zdrop
        attributes.heuristic.strategy = wf_heuristic_zdrop;
        attributes.heuristic.zdrop = MIN_OF_TWO(500, (int)(MIN_OF_TWO(plen, tlen) * 0.1));
        attributes.heuristic.steps_between_cutoffs = 100;
    }
    // Initialize Wavefront Aligner
    wavefront_aligner_t* const wf_aligner = wavefront_aligner_new(&attributes);
    uint8_t *p = pattern, *t = text;
    if (gap_aln == LONGCALLD_GAP_LEFT_ALN) { // reverse pattern and text
        p = (uint8_t*)malloc(plen); t = (uint8_t*)malloc(tlen);
        for (int i = 0; i < plen; ++i) p[i] = pattern[plen-i-1];
        for (int i = 0; i < tlen; ++i) t[i] = text[tlen-i-1];
    }

    // Align
    wfa_last_capped = 0;
    if (wfa_max_steps != INT_MAX) wavefront_aligner_set_max_alignment_steps(wf_aligner, wfa_max_steps);
    int wf_status = wavefront_align(wf_aligner, (const char*)p, plen, (const char*)t, tlen);
    if (wf_status != WF_STATUS_ALG_COMPLETED && wf_status != WF_STATUS_ALG_PARTIAL) {
        wfa_last_capped = (wf_status == WF_STATUS_MAX_STEPS_REACHED);
        if (cigar_length != NULL) *cigar_length = 0;
        if (cigar_buf != NULL) *cigar_buf = NULL;
        if (alg_length != NULL) *alg_length = 0;
        if (pattern_alg != NULL) *pattern_alg = NULL;
        if (text_alg != NULL) *text_alg = NULL;
        wavefront_aligner_delete(wf_aligner);
        if (gap_aln == LONGCALLD_GAP_LEFT_ALN) { free(p); free(t); }
        return -1;
    }
    cigar_t *cigar = wf_aligner->cigar;
    if (LONGCALLD_VERBOSE >= 3 && cigar != NULL && cigar->cigar_length > 0) {
        char *p_seq = (char*)malloc(plen+1); char *t_seq = (char*)malloc(tlen+1);
        for (int i = 0; i < plen; ++i) p_seq[i] = "ACGTN"[p[i]];
        for (int i = 0; i < tlen; ++i) t_seq[i] = "ACGTN"[t[i]];
        cigar_print_pretty(stderr, cigar, p_seq, plen, t_seq, tlen);
        free(p_seq); free(t_seq);
    }
    // collect cigar
    if (cigar_buf != NULL && cigar_length != NULL) {
        uint32_t *tmp_cigar;
        cigar_get_CIGAR(wf_aligner->cigar,true, &tmp_cigar, cigar_length);
        *cigar_buf = (uint32_t*)malloc(*cigar_length * sizeof(uint32_t));
        if (gap_aln == LONGCALLD_GAP_LEFT_ALN) { // reverse cigar
            for (int i = 0; i < *cigar_length; ++i) {
                (*cigar_buf)[i] = tmp_cigar[*cigar_length-i-1];
            }
        } else {
            for (int i = 0; i < *cigar_length; ++i) {
                (*cigar_buf)[i] = tmp_cigar[i];
            }
        }
    }
    // collect alignment string
    if (pattern_alg != NULL && text_alg != NULL) {
        *alg_length = wfa_collect_pretty_alignment(cigar, p, plen, t, tlen, pattern_alg, text_alg);
        if (gap_aln == LONGCALLD_GAP_LEFT_ALN) { // reverse alignment
            for (int i = 0; i < *alg_length/2; ++i) {
                uint8_t tmp = (*pattern_alg)[i];
                (*pattern_alg)[i] = (*pattern_alg)[*alg_length-i-1];
                (*pattern_alg)[*alg_length-i-1] = tmp;
                tmp = (*text_alg)[i];
                (*text_alg)[i] = (*text_alg)[*alg_length-i-1];
                (*text_alg)[*alg_length-i-1] = tmp;
            }
        }
    }
    // Free
    // fprintf(stderr, "%d vs %d Real time: %.3f sec. Score: %d\n", plen, tlen, realtime() - realtime0, cigar_score_gap_affine2p(cigar, &attributes.affine2p_penalties));
    wavefront_aligner_delete(wf_aligner); 
    if (gap_aln == LONGCALLD_GAP_LEFT_ALN) { free(p); free(t); }
    return 0;
}


int wfa_collect_diff_ins_seq(const call_var_opt_t *opt, uint8_t* large_seq, int large_len, uint8_t *small_seq, int small_len, uint8_t **diff_seq) {
    int diff_ins_len = 0;
    uint8_t *large_aln_seq = NULL, *small_aln_seq = NULL; int aln_len;
    wfa_end2end_aln(large_seq, large_len, small_seq, small_len,
                    opt->gap_aln, opt->mismatch, opt->gap_open1, opt->gap_ext1, opt->gap_open2, opt->gap_ext2, LONGCALLD_WFA_NO_HEURISTIC, LONGCALLD_WFA_AFFINE_2P, // no heuristic, affine-2p
                    NULL, NULL, &large_aln_seq, &small_aln_seq, &aln_len);
    // collect the largest insertion sequence
    int largest_ins_len = 0, largest_ins_pos = -1;
    for (int i = 0; i < aln_len; ++i) {
        if (small_aln_seq[i] == 5 && large_aln_seq[i] != 5) { // insertion in large_seq
            int ins_len = 0, j = i;
            while (j < aln_len && small_aln_seq[j] == 5 && large_aln_seq[j] != 5) {
                ins_len++; j++;
            }
            if (ins_len > largest_ins_len) {
                largest_ins_len = ins_len; largest_ins_pos = i;
            }
            i = j-1; // skip the inserted bases
        }
    }
    if (largest_ins_len > 0) {
        diff_ins_len = largest_ins_len;
        *diff_seq = (uint8_t*)malloc(largest_ins_len * sizeof(uint8_t));
        for (int i = 0; i < largest_ins_len; ++i) {
            (*diff_seq)[i] = large_aln_seq[largest_ins_pos+i];
        }
    }
    free(large_aln_seq);
    return diff_ins_len;
}

// trim query bases if longer than target
// update query_beg/query_end if query is shorter than target
void wfa_trim_aln_str(int full_cover, aln_str_t *aln_str) {
    if (LONGCALLD_NOISY_IS_NOT_COVER(full_cover) || LONGCALLD_NOISY_IS_BOTH_COVER(full_cover)) return;

    if ((LONGCALLD_NOISY_IS_LEFT_COVER(full_cover) && LONGCALLD_NOISY_IS_RIGHT_GAP(full_cover)) ||
        (LONGCALLD_NOISY_IS_RIGHT_COVER(full_cover) && LONGCALLD_NOISY_IS_LEFT_GAP(full_cover))) {
        aln_str->target_beg = 0; aln_str->target_end = aln_str->aln_len-1;
        aln_str->query_beg = 0; aln_str->query_end = aln_str->aln_len-1;
        return;
    }
    if (LONGCALLD_NOISY_IS_LEFT_COVER(full_cover)) { // left-cover: collect all aln_str until the very last = opertion
        int target_end = -1, query_end = -1;
        for (int i = aln_str->aln_len-1; i >= 0; --i) {
            if (query_end == -1) {
                if (aln_str->query_aln[i] != 5 && aln_str->target_aln[i] == aln_str->query_aln[i]) {
                    query_end = i;
                }
            }
            if (target_end == -1) {
                if (aln_str->target_aln[i] != 5) {
                    target_end = i;
                }
            }
            if (target_end != -1 && query_end != -1) break;
        }
        if (query_end == -1) query_end = target_end; // no = operation, query_end is the last position
        assert(query_end <= target_end);
        aln_str->aln_len = target_end+1;
        aln_str->target_beg = 0; aln_str->target_end = target_end;
        aln_str->query_beg = 0; aln_str->query_end = query_end;
        for (int i = query_end+1; i < aln_str->aln_len; ++i) {
            aln_str->query_aln[i] = 5;
        }
    } else if (LONGCALLD_NOISY_IS_RIGHT_COVER(full_cover)) { // right-cover: collect all aln_str until the very first = opertion
        int query_start = -1, target_start = -1;
        for (int i = 0; i < aln_str->aln_len; ++i) {
            if (query_start == -1) {
                if (aln_str->query_aln[i] != 5 && aln_str->target_aln[i] == aln_str->query_aln[i]) {
                    query_start = i;
                }
            }
            if (target_start == -1) {
                if (aln_str->target_aln[i] != 5) {
                    target_start = i;
                }
            }
            if (target_start != -1 && query_start != -1) break;
        }
        if (query_start == -1) query_start = target_start; // no = operation, query_start is the first position
        assert(query_start >= target_start);
        aln_str->aln_len = aln_str->aln_len - target_start;
        if (target_start != 0) {
            uint8_t *_target_aln = (uint8_t*)malloc(aln_str->aln_len * 2 * sizeof(uint8_t));
            for (int i = 0; i < aln_str->aln_len; ++i)
                _target_aln[i] = aln_str->target_aln[i+target_start];
            for (int i = 0; i < aln_str->aln_len; ++i)
                _target_aln[i+aln_str->aln_len] = aln_str->query_aln[i+target_start];
            free(aln_str->target_aln);
            aln_str->target_aln = _target_aln;
            aln_str->query_aln = _target_aln + aln_str->aln_len;
        }
        aln_str->target_beg = 0; aln_str->target_end = aln_str->aln_len-1;
        aln_str->query_beg = query_start-target_start; aln_str->query_end = aln_str->aln_len-1;
        for (int i = 0; i < aln_str->query_beg; ++i) {
            aln_str->query_aln[i] = 5;
        }
    }
}

// for read with full_cover as 1 or 2, collect beg/end positions of read mapped to target
int wfa_end2end_aln(uint8_t *pattern, int plen, uint8_t *text, int tlen,
                    int gap_aln, int b, int q, int e, int q2, int e2, int heuristic, int affine_gap,
                    uint32_t **cigar_buf, int *cigar_length, uint8_t **pattern_alg, uint8_t **text_alg, int *alg_length) {
    return wfa_end2end_aln_mem(0, 0, pattern, plen, text, tlen, gap_aln, b, q, e, q2, e2, heuristic, affine_gap,
                               cigar_buf, cigar_length, pattern_alg, text_alg, alg_length);
}

int wfa_collect_aln_str(const call_var_opt_t *opt, uint8_t *target, int tlen, uint8_t *query, int qlen, int full_cover, int heuristic, int affine_gap, aln_str_t *aln_str) {
    if (LONGCALLD_NOISY_IS_NOT_COVER(full_cover)) return 0;
    aln_str->target_aln = 0; aln_str->query_aln = 0; aln_str->aln_len = 0;
    int gap_aln = opt->gap_aln, b = opt->mismatch, q = opt->gap_open1, e = opt->gap_ext1, q2 = opt->gap_open2, e2 = opt->gap_ext2;
    if (LONGCALLD_NOISY_IS_BOTH_COVER(full_cover)) {
        wfa_end2end_aln(target, tlen, query, qlen, gap_aln, b, q, e, q2, e2, heuristic, affine_gap,
                        NULL, NULL, &aln_str->target_aln, &aln_str->query_aln, &aln_str->aln_len);
        aln_str->target_beg = 0; aln_str->target_end = aln_str->aln_len-1;
        aln_str->query_beg = 0; aln_str->query_end = aln_str->aln_len-1;
    } else { // trim aln_str if query is longer than target; 
             // update query_beg/query_end if query is shorter than target
        int ext_direction = (LONGCALLD_NOISY_IS_LEFT_COVER(full_cover)) ? LONGCALLD_EXT_ALN_LEFT_TO_RIGHT : LONGCALLD_EXT_ALN_RIGHT_TO_LEFT;
        // trim target/query to opt->LONGCALLD_PARTIAL_ALN_RATIO
        double ratio = opt->partial_aln_ratio;
        int _tlen = tlen, _qlen = qlen;
        int _t_start = 0, _q_start = 0;
        if (ext_direction == LONGCALLD_EXT_ALN_LEFT_TO_RIGHT) {
            if (tlen > qlen*ratio) _tlen = qlen*ratio;
            else if (qlen > tlen*ratio) _qlen = tlen*ratio;
        } else { // Right to Left
            if (tlen > qlen*ratio) { _tlen = qlen*ratio; _t_start = tlen - _tlen; }
            else if (qlen > tlen*ratio) { _qlen = tlen*ratio; _q_start = qlen - _qlen; }
        }
        if (ext_direction == LONGCALLD_EXT_ALN_LEFT_TO_RIGHT) {
            gap_aln = (gap_aln == LONGCALLD_GAP_LEFT_ALN) ? LONGCALLD_GAP_RIGHT_ALN : LONGCALLD_GAP_LEFT_ALN;
        }
        // partial alignment
        wfa_end2end_aln(target+_t_start, _tlen, query+_q_start, _qlen, 
                        gap_aln, b, q, e, q2, e2, LONGCALLD_WFA_ZDROP, LONGCALLD_WFA_AFFINE_2P,
                        NULL, NULL, &aln_str->target_aln, &aln_str->query_aln, &aln_str->aln_len);
        // do not trim for end-gaps
        wfa_trim_aln_str(full_cover, aln_str);
        // target_beg/target_end/query_beg/query_end are alignment columns and target positions are counted from
        // column 0, so the skipped target prefix is put back as gap-only columns (the skipped query prefix needs none)
        if (_t_start > 0 && aln_str->target_aln != NULL && aln_str->aln_len > 0) {
            int new_len = aln_str->aln_len + _t_start;
            uint8_t *_target_aln = (uint8_t*)malloc(new_len * 2 * sizeof(uint8_t));
            for (int i = 0; i < _t_start; ++i) { _target_aln[i] = target[i]; _target_aln[new_len+i] = 5; }
            for (int i = 0; i < aln_str->aln_len; ++i) {
                _target_aln[_t_start+i] = aln_str->target_aln[i];
                _target_aln[new_len+_t_start+i] = aln_str->query_aln[i];
            }
            free(aln_str->target_aln);
            aln_str->target_aln = _target_aln; aln_str->query_aln = _target_aln + new_len; aln_str->aln_len = new_len;
            aln_str->target_beg = 0; aln_str->target_end = new_len-1;
            aln_str->query_beg += _t_start; aln_str->query_end += _t_start;
        }
    }
    if (LONGCALLD_VERBOSE >= 2) {
        fprintf(stderr, ">target %d:%d-%d\n", tlen, aln_str->target_beg, aln_str->target_end);
        for (int i = 0; i < aln_str->aln_len; ++i) fprintf(stderr, "%c", "ACGTN-"[aln_str->target_aln[i]]);
        fprintf(stderr, "\n>query %d:%d-%d\n", qlen, aln_str->query_beg, aln_str->query_end);
        for (int i = 0; i < aln_str->aln_len; ++i) fprintf(stderr, "%c", "ACGTN-"[aln_str->query_aln[i]]);
        fprintf(stderr, "\n");
    }
    return 0;
}

int end2end_aln(const call_var_opt_t *opt, char *tseq, int tlen, uint8_t *qseq, int qlen, uint32_t **cigar_buf) {
    if (qlen <= 0 || tlen <= 0) return 0;
    // int min_len = MIN_OF_TWO(tlen, qlen)
    int max_len = MAX_OF_TWO(tlen, qlen);
    // int delta_len = MAX_OF_TWO(1, max_len - min_len);

    uint8_t *tseq2 = (uint8_t*)malloc(max_len);
    for (int i = 0; i < tlen; ++i) tseq2[i] = nst_nt4_table[(uint8_t)tseq[i]];
    // use wfa if the length difference is small
    // if (max_len * delta_len + delta_len * delta_len < max_len * min_len) {
    int cigar_len = 0;
    wfa_end2end_aln(tseq2, tlen, qseq, qlen, opt->gap_aln, opt->mismatch, opt->gap_open1, opt->gap_ext1, opt->gap_open2, opt->gap_ext2, LONGCALLD_WFA_NO_HEURISTIC, LONGCALLD_WFA_AFFINE_2P,
                    cigar_buf, &cigar_len, NULL, NULL, NULL);
    // } else { // if (max_len - min_len > 1000) { // use ksw2 if the length difference is large
        // cigar_len = ksw2_aln(opt->gap_aln, tseq2, tlen, qseq, qlen, opt->match, opt->mismatch, opt->gap_open1, opt->gap_ext1, opt->gap_open2, opt->gap_ext2, cigar_buf);
    // }
    free(tseq2);
    return cigar_len;
}

void collect_aln_beg_end(uint32_t *cigar_buf, int cigar_len, int ext_direction, int ref_len, int *ref_beg, int *ref_end, int read_len, int *read_beg, int *read_end) {
    *ref_beg = 1; *read_beg = 1; *ref_end = ref_len, *read_end = read_len;
    if (ext_direction == LONGCALLD_EXT_ALN_LEFT_TO_RIGHT) {
        int tmp_ref_end = 0, tmp_read_end = 0;
        for (int i = 0; i < cigar_len; ++i) { // collect ref_end & read_end with last M cigar
            int op = cigar_buf[i] & 0xf; int len = cigar_buf[i] >> 4;
            if (op == BAM_CEQUAL || op == BAM_CMATCH) {
                tmp_ref_end += len; tmp_read_end += len;
                *ref_end = tmp_ref_end; *read_end = tmp_read_end;
            } else if (op == BAM_CDIFF) {
                tmp_ref_end += len; tmp_read_end += len;
            } else if (op == BAM_CDEL) {
                tmp_ref_end += len;
            } else if (op == BAM_CINS) {
                tmp_read_end += len;
            } else continue;
        }
    } else {
        int tmp_ref_beg = ref_len+1, tmp_read_beg = read_len+1;
        for (int i=cigar_len-1; i>=0; i--) {
            int op = cigar_buf[i] & 0xf; int len = cigar_buf[i] >> 4;
            if (op == BAM_CEQUAL || op == BAM_CMATCH) {
                tmp_ref_beg -= len; tmp_read_beg -= len;
                *ref_beg = tmp_ref_beg; *read_beg = tmp_read_beg;
            } else if (op == BAM_CDIFF) {
                tmp_ref_beg -= len; tmp_read_beg -= len;
            } else if (op == BAM_CDEL) {
                tmp_ref_beg -= len;
            } else if (op == BAM_CINS) {
                tmp_read_beg -= len;
            } else continue;
        }
    }
}

// calculate the beg/end positions of query mapped to target using WFA
// limit to (partial_aln_ratio*short_len vs short_len)
int cal_wfa_partial_aln_beg_end(int ext_direction, const call_var_opt_t *opt, uint8_t *_target, int _tlen, uint8_t *_query, int _qlen, int *target_beg, int *target_end, int *query_beg, int *query_end) {
    int gap_aln = opt->gap_aln, b = opt->mismatch, q = opt->gap_open1, e = opt->gap_ext1, q2 = opt->gap_open2, e2 = opt->gap_ext2;
    double ratio = opt->partial_aln_ratio;
    int tlen = _tlen, qlen = _qlen;
    uint8_t *target = _target, *query = _query;
    if (ext_direction == LONGCALLD_EXT_ALN_LEFT_TO_RIGHT) {
        if (_tlen > _qlen * ratio) tlen = (int)(_qlen * ratio);
        else if (_qlen > _tlen * ratio) qlen = (int)(_tlen * ratio);
    } else if (ext_direction == LONGCALLD_EXT_ALN_RIGHT_TO_LEFT) {
        if (_tlen > _qlen * ratio) {
            target = _target + _tlen - (int)(_qlen * ratio);
            tlen = (int)(_qlen * ratio);
        } else if (_qlen > _tlen * ratio) {
            query = _query + _qlen - (int)(_tlen * ratio);
            qlen = (int)(_tlen * ratio);
        }
    }
    uint32_t *cigar_buf = NULL;
    int ret = 1, cigar_len;
    // left to right
    if (ext_direction == LONGCALLD_EXT_ALN_LEFT_TO_RIGHT) {
        gap_aln = (gap_aln == LONGCALLD_GAP_RIGHT_ALN) ? LONGCALLD_GAP_LEFT_ALN : LONGCALLD_GAP_RIGHT_ALN;
    }
    // fprintf(stderr, "WFA partial aln: %d vs %d (Raw: %d vs %d)\n", tlen, qlen, _tlen, _qlen);

    int min_len = MIN_OF_TWO(tlen, qlen);
    if (ext_direction == LONGCALLD_EXT_ALN_LEFT_TO_RIGHT) {
        int x_gaps = edlib_xgaps(target, min_len, query, min_len);
        if (x_gaps > min_len * 0.10) { return 0;
        }
    } else {
        int x_gaps = edlib_xgaps(target + tlen - min_len, min_len, query + qlen - min_len, min_len);
        if (x_gaps > min_len * 0.10) { return 0;
        }
    }
    wfa_end2end_aln(target, tlen, query, qlen, gap_aln, b, q, e, q2, e2, LONGCALLD_WFA_NO_HEURISTIC, LONGCALLD_WFA_AFFINE_2P, &cigar_buf, &cigar_len, NULL, NULL, NULL);
    if (cigar_len == 0) ret = 0;
    else collect_aln_beg_end(cigar_buf, cigar_len, ext_direction, _tlen, target_beg, target_end, _qlen, query_beg, query_end);
    if (cigar_buf != NULL) free(cigar_buf);
    return ret;
}

int collect_partial_aln_beg_end(const call_var_opt_t *opt, int sampling_reads,
                                uint8_t *target, int tlen, int target_full_cover, uint8_t *query, int qlen, int query_full_cover, 
                                int *target_beg, int *target_end, int *query_beg, int *query_end) {
    *target_beg = 1, *target_end = tlen, *query_beg = 1, *query_end = qlen;
    int ret = 1;
    assert(LONGCALLD_NOISY_IS_BOTH_COVER(target_full_cover) != 0);
    if (LONGCALLD_NOISY_IS_BOTH_COVER(target_full_cover)) { // ref is full cover
        if (LONGCALLD_NOISY_IS_BOTH_COVER(query_full_cover) || 
           (LONGCALLD_NOISY_IS_LEFT_COVER(query_full_cover) && LONGCALLD_NOISY_IS_RIGHT_GAP(query_full_cover)) ||
           (LONGCALLD_NOISY_IS_RIGHT_COVER(query_full_cover) && LONGCALLD_NOISY_IS_LEFT_GAP(query_full_cover))) {
            if (sampling_reads) { // skip full-cover reads with >5% ED in sampling mode
                // int edit_dis = edlib_edit_distance(target, tlen, query, qlen);
                // if (edit_dis > MIN_OF_TWO(tlen, qlen) * 0.05) {
                int x_gaps = edlib_xgaps(target, tlen, query, qlen);
                if (x_gaps > MIN_OF_TWO(tlen, qlen) * 0.10) {
                    // if (LONGCALLD_VERBOSE >= 3) fprintf(stderr, "Skipped in POA due to high edit distance (%d vs %d): %d\n", tlen, qlen, edit_dis);
                    if (LONGCALLD_VERBOSE >= 3) fprintf(stderr, "Skipped in POA due to high # of X-gaps (%d vs %d): %d\n", tlen, qlen, x_gaps);
                    return 0;
                }
            } 
            return 1;
        } else {
            if (LONGCALLD_NOISY_IS_LEFT_COVER(query_full_cover)) {
                ret = cal_wfa_partial_aln_beg_end(LONGCALLD_EXT_ALN_LEFT_TO_RIGHT, opt, target, tlen, query, qlen, target_beg, target_end, query_beg, query_end);
            } else if (LONGCALLD_NOISY_IS_RIGHT_COVER(query_full_cover)) {
                ret = cal_wfa_partial_aln_beg_end(LONGCALLD_EXT_ALN_RIGHT_TO_LEFT, opt, target, tlen, query, qlen, target_beg, target_end, query_beg, query_end);
            }
        }
    } else if (LONGCALLD_NOISY_IS_LEFT_COVER(target_full_cover)) { // ref is left-cover
        if (LONGCALLD_NOISY_IS_LEFT_COVER(query_full_cover) == 0) _err_error_exit("Target is left-cover but read is not left-cover\n");
        ret = cal_wfa_partial_aln_beg_end(LONGCALLD_EXT_ALN_LEFT_TO_RIGHT, opt, target, tlen, query, qlen, target_beg, target_end, query_beg, query_end);
    } else if (LONGCALLD_NOISY_IS_RIGHT_COVER(target_full_cover)) { // ref is right-cover
        if (LONGCALLD_NOISY_IS_RIGHT_COVER(query_full_cover) == 0) _err_error_exit("Ref is right-cover but read is not right-cover\n");
        ret = cal_wfa_partial_aln_beg_end(LONGCALLD_EXT_ALN_RIGHT_TO_LEFT, opt, target, tlen, query, qlen, target_beg, target_end, query_beg, query_end);
    } else _err_error_exit("Target is not left or right-cover\n");
    return ret;
}

int collect_exc_beg_end(const call_var_opt_t *opt, abpoa_t *ab, abpoa_para_t *abpt, int sampling_reads, int n_reads, char **names, uint8_t **read_seqs, 
                        int *read_lens, int *read_full_cover, int *exc_begs, int *exc_ends, int *seq_beg_cuts, int *seq_end_cuts) {
    for (int read_i = 1; read_i < n_reads; ++read_i) {
        int ref_beg, ref_end, read_beg, read_end, beg_id, end_id;
        if (collect_partial_aln_beg_end(opt, sampling_reads, read_seqs[0], read_lens[0], read_full_cover[0], read_seqs[read_i], read_lens[read_i], read_full_cover[read_i], &ref_beg, &ref_end, &read_beg, &read_end) == 0) {
            exc_begs[read_i] = -1, exc_ends[read_i] = -1;
            continue;
        }
        beg_id = ref_beg+1, end_id = ref_end+1;
        seq_beg_cuts[read_i] = read_beg - 1, seq_end_cuts[read_i] = read_lens[read_i] - read_end;
        abpoa_subgraph_nodes(ab, abpt, beg_id, end_id, exc_begs+read_i, exc_ends+read_i);
    }
    return 0;
}

int abpoa_partial_aln_msa_cons(const call_var_opt_t *opt, abpoa_t *ab, int sampling_reads, int n_reads, int *read_ids, uint8_t **read_seqs, uint8_t **read_quals, int *read_lens, int *read_full_cover, char **names,
                               int max_n_cons, int *cons_lens, uint8_t **cons_seqs, int *clu_n_seqs, int **clu_read_ids, int *msa_seq_lens, uint8_t **msa_seqs) {
    // abpoa_t *ab = abpoa_init();
    int needs_free_ab = 0;
    if (ab == NULL) {
        ab = abpoa_init(); needs_free_ab = 1;
    }
    abpoa_para_t *abpt = abpoa_init_para();
    // if (opt->out_somatic) abpt->wf = 0.01; // more accurate for somatic variant calling
    // else abpt->wf = 0.001; // limit memory usage for long sequences
    abpt->cons_algrm = ABPOA_MF;
    abpt->sub_aln = 1;
    abpt->inc_path_score = 1;
    if (cons_lens != NULL && cons_seqs != NULL) abpt->out_cons = 1; else abpt->out_cons = 0;
    if (msa_seq_lens != NULL && msa_seqs != NULL) abpt->out_msa = 1; else abpt->out_msa = 0;
    abpt->max_n_cons = max_n_cons; abpt->min_freq = opt->min_af;
    abpt->match = opt->match; abpt->mismatch = opt->mismatch;
    abpt->gap_open1 = opt->gap_open1; abpt->gap_ext1 = opt->gap_ext1;
    abpt->gap_open2 = opt->gap_open2; abpt->gap_ext2 = opt->gap_ext2;

    if (LONGCALLD_VERBOSE >= 2) abpt->out_msa = 1;
    abpoa_post_set_para(abpt);
    ab->abs->n_seq = n_reads;
    for (int i = 0; i < n_reads; ++i) {
        if (LONGCALLD_VERBOSE >= 3) {
            fprintf(stderr, ">%s %d %d\n", names[i], read_lens[i], read_full_cover[i]);
            for (int j = 0; j < read_lens[i]; ++j) {
                fprintf(stderr, "%c", "ACGTN"[read_seqs[i][j]]);
            } fprintf(stderr, "\n");
        }
        abpoa_res_t res; res.graph_cigar = 0, res.n_cigar = 0;
        int exc_beg = 0, exc_end = 1, seq_beg_cut = 0, seq_end_cut = 0;
        if (i != 0) {
            int ref_beg, ref_end, read_beg, read_end, beg_id, end_id;
            if (collect_partial_aln_beg_end(opt, sampling_reads, read_seqs[0], read_lens[0], read_full_cover[0], read_seqs[i], read_lens[i], read_full_cover[i], &ref_beg, &ref_end, &read_beg, &read_end) == 0) continue;
            beg_id = ref_beg+1, end_id = ref_end+1;
            seq_beg_cut = read_beg - 1, seq_end_cut = read_lens[i] - read_end;
            abpoa_subgraph_nodes(ab, abpt, beg_id, end_id, &exc_beg, &exc_end);
            // fprintf(stderr, "beg_id: %d, end_id: %d, read_beg: %d, read_end: %d, exc_beg: %d, exc_end: %d\n", beg_id, end_id, read_beg, read_end, exc_beg, exc_end);
        }
        if (LONGCALLD_VERBOSE >= 3) fprintf(stderr, "ExcBeg: %d, ExcEnd: %d, SeqBegCut: %d, SeqEndCut: %d, FullCover: %d\n", exc_beg, exc_end, seq_beg_cut, seq_end_cut, read_full_cover[i]);
        abpoa_align_sequence_to_subgraph(ab, abpt, exc_beg, exc_end, read_seqs[i]+seq_beg_cut, read_lens[i]-seq_beg_cut-seq_end_cut, &res);
        // abpoa_add_subgraph_alignment(ab, abpt, exc_begs[i], exc_ends[i], read_seqs[i]+seq_beg_cuts[i], NULL, read_lens[i]-seq_beg_cuts[i]-seq_end_cuts[i], NULL, res, i, n_reads, 1);
        abpoa_add_subgraph_alignment(ab, abpt, exc_beg, exc_end, read_seqs[i]+seq_beg_cut, NULL, read_lens[i]-seq_beg_cut-seq_end_cut, NULL, res, i, n_reads, 0);
        if (res.n_cigar) free(res.graph_cigar);
    }
    if (LONGCALLD_VERBOSE >= 2) abpoa_output(ab, abpt, stderr);
    else abpoa_output(ab, abpt, NULL);
    abpoa_cons_t *abc = ab->abc;
    
    int n_cons = 0;
    // cons bases
    if (cons_lens != NULL && cons_seqs != NULL) {
        if (abc->n_cons > 0) {
            for (int i = 0; i < abc->n_cons; ++i) {
                cons_lens[i] = abc->cons_len[i];
                cons_seqs[i] = (uint8_t*)malloc(abc->cons_len[i] * sizeof(uint8_t));
                for (int j = 0; j < abc->cons_len[i]; ++j) {
                    cons_seqs[i][j] = abc->cons_base[i][j];
                }
                if (LONGCALLD_VERBOSE >= 2) fprintf(stderr, "ConsLen: %d\n", abc->cons_len[i]);
            }
            if (clu_n_seqs != NULL && clu_read_ids != NULL) {
                if (abc->n_cons > 1) {
                    for (int i = 0; i < abc->n_cons; ++i) {
                        clu_n_seqs[i] = abc->clu_n_seq[i];
                        clu_read_ids[i] = (int*)malloc(abc->clu_n_seq[i] * sizeof(int));
                        for (int j = 0; j < abc->clu_n_seq[i]; ++j) {
                            clu_read_ids[i][j] = read_ids[abc->clu_read_ids[i][j]];
                        }
                    }
                } else {
                    *clu_n_seqs = n_reads;
                    *clu_read_ids = (int*)malloc(n_reads * sizeof(int));
                    for (int i = 0; i < n_reads; ++i) (*clu_read_ids)[i] = read_ids[i];
                }
            }
            n_cons = abc->n_cons;
        } else {
            if (LONGCALLD_VERBOSE >= 2) fprintf(stderr, "Unable to call consensus: %s\n", names[0]);
        }
    }

    // msa bases, include consensus sequences
    if (msa_seq_lens != NULL && msa_seqs != NULL) {
        *msa_seq_lens = abc->msa_len;
        for (int i = 0; i < abc->n_seq+n_cons; ++i) {
            msa_seqs[i] = (uint8_t*)malloc(abc->msa_len * sizeof(uint8_t));
            for (int j = 0; j < abc->msa_len; ++j) {
                msa_seqs[i][j] = abc->msa_base[i][j];
            }
        }
    }
    abpoa_free_para(abpt); if (needs_free_ab) abpoa_free(ab);
    return n_cons;
 }

// TODO: for large deletions
// int two_sides_abpoa_partial_aln_msa_cons(const call_var_opt_t *opt, abpoa_t *ab, int n_reads, int *read_ids, uint8_t **read_seqs, uint8_t **read_quals, int *read_lens, int *read_full_cover, char **names,
//                                          int max_n_cons, int *cons_lens, uint8_t **cons_seqs, int *clu_n_seqs, int **clu_read_ids) {
//     int n_cons = 0;
//     // all left-cover: manually set full-cover of first read as 3, collect consensus
//     int n_left_cover = 0;
//     // all right-cover: manually set full-cover of first read as 3, collect consensus
//     int n_right_cover = 0;
//     // merge consensus
//     return n_cons;
//  }

  // XXX limit abpoa memory usage, avoid memory allocation failure
 int abpoa_aln_msa_cons(const call_var_opt_t *opt, int n_reads, int *read_ids, uint8_t **read_seqs, int *read_lens, int max_n_cons,
                        int *cons_lens, uint8_t **cons_seqs,
                        int *clu_n_seqs, int **clu_read_ids, int *msa_seq_len, uint8_t ***msa_seq) {
    abpoa_t *ab = abpoa_init();
    abpoa_para_t *abpt = abpoa_init_para();
    // if (opt->out_somatic) abpt->wf = 0.01; // XXX for more accurate somatic variant calling
    // else abpt->wf = 0.001; // limit memory usage for long sequences
    abpt->wb = -1;
    abpt->inc_path_score = 1;
    abpt->out_msa = 0; abpt->out_cons = 1;
    if (msa_seq != NULL && msa_seq_len != NULL) abpt->out_msa = 1; else abpt->out_msa = 0;
    if (LONGCALLD_VERBOSE >= 2) abpt->out_msa = 1;
    abpt->cons_algrm = ABPOA_MF;
    abpt->max_n_cons = max_n_cons; abpt->min_freq = opt->min_af;
    abpt->match = opt->match; abpt->mismatch = opt->mismatch;
    abpt->gap_open1 = opt->gap_open1; abpt->gap_ext1 = opt->gap_ext1;
    abpt->gap_open2 = opt->gap_open2; abpt->gap_ext2 = opt->gap_ext2;
    abpoa_post_set_para(abpt);

    if (LONGCALLD_VERBOSE >= 2) {
        fprintf(stderr, "For abPOA (max %d cons, min_freq: %.2f): %d\n", max_n_cons, abpt->min_freq, n_reads);
        abpoa_msa(ab, abpt, n_reads, NULL, read_lens, read_seqs, NULL, stderr);
    } else abpoa_msa(ab, abpt, n_reads, NULL, read_lens, read_seqs, NULL, NULL);
    abpoa_cons_t *abc = ab->abc;
    
    int n_cons = 0;
    // cons bases
    if (abc->n_cons > 0) {
        for (int i = 0; i < abc->n_cons; ++i) {
            cons_lens[i] = abc->cons_len[i];
            cons_seqs[i] = (uint8_t*)malloc(abc->cons_len[i] * sizeof(uint8_t));
            for (int j = 0; j < abc->cons_len[i]; ++j) {
                cons_seqs[i][j] = abc->cons_base[i][j];
            }
            if (LONGCALLD_VERBOSE >= 2) fprintf(stderr, "ConsLen: %d\n", abc->cons_len[i]);
        }
        if (abc->n_cons == 2) { // collect read ids for each cluster
            for (int i = 0; i < 2; ++i) {
                clu_n_seqs[i] = abc->clu_n_seq[i];
                clu_read_ids[i] = (int*)malloc(abc->clu_n_seq[i] * sizeof(int));
                for (int j = 0; j < abc->clu_n_seq[i]; ++j) {
                    clu_read_ids[i][j] = read_ids[abc->clu_read_ids[i][j]];
                }
            }
        } else {
            *clu_n_seqs = n_reads;
            *clu_read_ids = (int*)malloc(n_reads * sizeof(int));
            for (int i = 0; i < n_reads; ++i) (*clu_read_ids)[i] = read_ids[i];
        }
        n_cons = abc->n_cons;
    } else {
        _err_error_exit("No Consensus.\n");
    }
    // msa bases
    if (msa_seq != NULL && msa_seq_len != NULL) {
        // split msa and cons into 2 clusers
        size_t msize = abc->msa_len * sizeof(uint8_t);
        for (int i = 0; i < n_cons; ++i) {
            msa_seq_len[i] = abc->msa_len;
            for (int j = 0; j < abc->clu_n_seq[i]; ++j) {
                msa_seq[i][j] = (uint8_t*)malloc(msize);
                for (int k = 0; k < abc->msa_len; ++k)
                    msa_seq[i][j][k] = abc->msa_base[abc->clu_read_ids[i][j]][k];
            }
            msa_seq[i][abc->clu_n_seq[i]] = (uint8_t*)malloc(msize);
            for (int j = 0; j < abc->msa_len; ++j)
                msa_seq[i][abc->clu_n_seq[i]][j] = abc->msa_base[abc->n_seq+i][j];
        }
    }
    abpoa_free_para(abpt); abpoa_free(ab); 
    return n_cons;
}

int full_cover_cmp(int cover1, int cover2) {
    if (cover1 == cover2) return 0;
    if (LONGCALLD_NOISY_IS_BOTH_COVER(cover1)) return 1;
    else if (LONGCALLD_NOISY_IS_BOTH_COVER(cover2)) return -1;
    if (LONGCALLD_NOISY_IS_LEFT_COVER(cover1) && LONGCALLD_NOISY_IS_LEFT_COVER(cover2)) return 0;
    if (LONGCALLD_NOISY_IS_RIGHT_COVER(cover1) && LONGCALLD_NOISY_IS_RIGHT_COVER(cover2)) return 0;
    return cover1-cover2;
}

// sort reads in noisy regions based on full-cover, [error-rate], read-length; error-rate is optional
void sort_noisy_region_reads(int n_reads, int *read_ids, int *read_lens, uint8_t **read_seqs, uint8_t **read_quals, uint8_t *strands, int *full_covers, char **names,
                             int *haps, hts_pos_t *phase_sets, int use_error_rate) {
    double *read_error_rates = NULL;
    if (use_error_rate) {
        read_error_rates = (double*)malloc(n_reads * sizeof(double));
        for (int i = 0; i < n_reads; ++i)
            read_error_rates[i] = calc_read_error_rate(read_lens[i], read_quals[i]);
    }
    for (int i = 0; i < n_reads-1; ++i) {
        for (int j = i+1; j < n_reads; ++j) {
            int cover_cmp = full_cover_cmp(full_covers[i], full_covers[j]);
            if (cover_cmp < 0 || 
               (cover_cmp == 0 && use_error_rate && read_error_rates[i] > read_error_rates[j]) ||
               (cover_cmp == 0 && ((use_error_rate && read_error_rates[i] == read_error_rates[j]) || !use_error_rate) && read_lens[i] < read_lens[j])) {
            // if (full_covers[i] < full_covers[j] || (full_covers[i] == full_covers[j] && read_lens[i] < read_lens[j])) {
                int tmp_full_cover = full_covers[i]; full_covers[i] = full_covers[j]; full_covers[j] = tmp_full_cover;
                int tmp_id = read_ids[i]; read_ids[i] = read_ids[j]; read_ids[j] = tmp_id;
                int tmp_len = read_lens[i]; read_lens[i] = read_lens[j]; read_lens[j] = tmp_len;
                uint8_t *tmp_seq = read_seqs[i]; read_seqs[i] = read_seqs[j]; read_seqs[j] = tmp_seq;
                uint8_t *tmp_qual = read_quals[i]; read_quals[i] = read_quals[j]; read_quals[j] = tmp_qual;
                uint8_t tmp_strand = strands[i]; strands[i] = strands[j]; strands[j] = tmp_strand;
                char *tmp_name = names[i]; names[i] = names[j]; names[j] = tmp_name;
                int tmp_hap = haps[i]; haps[i] = haps[j]; haps[j] = tmp_hap;
                hts_pos_t tmp_ps = phase_sets[i]; phase_sets[i] = phase_sets[j]; phase_sets[j] = tmp_ps;
                double tmp_error_rate = 0.0;
                if (use_error_rate) {
                    tmp_error_rate = read_error_rates[i]; read_error_rates[i] = read_error_rates[j]; read_error_rates[j] = tmp_error_rate;
                }
            }
        }
    }
    if (use_error_rate) free(read_error_rates);
}

static int add_phase_set(hts_pos_t ps, hts_pos_t *uniq_phase_sets, int *n_uniq_phase_sets) {
    int i;
    for (i = 0; i < *n_uniq_phase_sets; ++i) {
        if (uniq_phase_sets[i] == ps) return i;
    }
    uniq_phase_sets[i] = ps;
    (*n_uniq_phase_sets)++;
    return i;
}

// hp_start/end: flank_base
int is_homopolymer(uint8_t *seq, int seq_len, int flank_len, int *hp_start, int *hp_end, int *hp_len) {
    if (seq_len < 2*flank_len || seq_len > 2*flank_len + 50) return 0; // skip too short or too long cons
    int min_hp_len = 5; *hp_start = -1, *hp_end = -1, *hp_len = 0;
    for (int i = flank_len-1; i < seq_len-flank_len+1; ++i) {
        if (seq[i] == seq[i-1]) {
            if (*hp_start == -1) *hp_start = i-2;
            (*hp_len)++;
        } else {
            if (*hp_len >= min_hp_len) {
                *hp_end = i;
                return 1;
            } else {
                *hp_start = -1; *hp_len = 0;
            }
        }
    }
    if (*hp_len >= min_hp_len) {
        *hp_end = *hp_start + *hp_len + 1;
        return 1;
    }
    return 0;
}

// construct cons vs read alignment string based on abPOA's msa result
// read: A-CCGGT--A
// cons: ACCC-GT--A
// output:
// read: A-CCGGTA
// cons: ACCC-GTA
int make_cons_read_aln_str(const call_var_opt_t *opt, uint8_t *cons_str, uint8_t *read_str, int msa_len, int full_cover, aln_str_t *cons_read_aln_str) {
    int aln_len = 0;
    cons_read_aln_str->aln_len = 0;
    cons_read_aln_str->target_aln = (uint8_t*)malloc(msa_len * 2 * sizeof(uint8_t));
    cons_read_aln_str->query_aln = cons_read_aln_str->target_aln + msa_len;
    for (int i = 0; i < msa_len; ++i) {
        if (read_str[i] != 5 || cons_str[i] != 5) { // not both gaps
            cons_read_aln_str->target_aln[aln_len] = cons_str[i];
            cons_read_aln_str->query_aln[aln_len] = read_str[i];
            aln_len++;
        }
    }
    cons_read_aln_str->aln_len = aln_len;
    cons_read_aln_str->target_beg = 0; cons_read_aln_str->target_end = aln_len - 1;
    cons_read_aln_str->query_beg = 0; cons_read_aln_str->query_end = aln_len - 1;
    wfa_trim_aln_str(full_cover, cons_read_aln_str);

    if (LONGCALLD_VERBOSE >= 2) {
        fprintf(stderr, ">target %d-%d\n", cons_read_aln_str->target_beg, cons_read_aln_str->target_end);
        for (int i = 0; i < cons_read_aln_str->aln_len; ++i) fprintf(stderr, "%c", "ACGTN-"[cons_read_aln_str->target_aln[i]]);
        fprintf(stderr, "\n>query %d-%d\n", cons_read_aln_str->query_beg, cons_read_aln_str->query_end);
        for (int i = 0; i < cons_read_aln_str->aln_len; ++i) fprintf(stderr, "%c", "ACGTN-"[cons_read_aln_str->query_aln[i]]);
        fprintf(stderr, "\n");
    }
    return aln_len;
}

int make_ref_read_aln_str(const call_var_opt_t *opt, aln_str_t *ref_cons_aln_str, aln_str_t *cons_read_aln_str, aln_str_t *ref_read_aln_str) {
    int aln_len = 0;
    int max_msa_len = ref_cons_aln_str->aln_len + cons_read_aln_str->aln_len;
    ref_read_aln_str->aln_len = 0;
    ref_read_aln_str->target_aln = (uint8_t*)malloc(max_msa_len * 2 * sizeof(uint8_t));
    ref_read_aln_str->query_aln = ref_read_aln_str->target_aln + max_msa_len;
    int i = 0, j = 0;
    while (i < ref_cons_aln_str->aln_len && j < cons_read_aln_str->aln_len) {
        if (ref_cons_aln_str->query_aln[i] == 5 && cons_read_aln_str->target_aln[j] == 5) { // both are gaps: extract the seqs for ref and read, until the non-gap base of cons, do WFA
            int ref_del_len = 1, read_del_len = 1;
            while (i+ref_del_len < ref_cons_aln_str->aln_len && ref_cons_aln_str->query_aln[i+ref_del_len] == 5) ref_del_len++;
            while (j+read_del_len < cons_read_aln_str->aln_len && cons_read_aln_str->target_aln[j+read_del_len] == 5) read_del_len++;
            uint8_t *ref_aln=0, *read_aln=0; int del_aln_len;
            // print ref_cons_aln_str->target_aln[i] and cons_read_aln_str->query_aln[j] for debugging
            if (LONGCALLD_VERBOSE >= 3) {
                fprintf(stderr, "Ref-Cons Gap: %d, Cons-Read Gap: %d\n", ref_del_len, read_del_len);
                fprintf(stderr, "Ref-Cons Target Aln: ");
                for (int k = i; k < i+ref_del_len; ++k) {
                    fprintf(stderr, "%c", "ACGTN-"[ref_cons_aln_str->target_aln[k]]);
                } fprintf(stderr, "\n");
                fprintf(stderr, "Cons-Read Query Aln: ");
                for (int k = j; k < j+read_del_len; ++k) {
                    fprintf(stderr, "%c", "ACGTN-"[cons_read_aln_str->query_aln[k]]);
                } fprintf(stderr, "\n");
            }
            wfa_end2end_aln(ref_cons_aln_str->target_aln+i, ref_del_len, cons_read_aln_str->query_aln+j, read_del_len, 
                            opt->gap_aln, opt->mismatch, opt->gap_open1, opt->gap_ext1, opt->gap_open2, opt->gap_ext2, LONGCALLD_WFA_NO_HEURISTIC, LONGCALLD_WFA_AFFINE_2P, // no heuristic, affine-2p
                            NULL, NULL, &ref_aln, &read_aln, &del_aln_len);
            for (int k = 0; k < del_aln_len; ++k) {
                ref_read_aln_str->target_aln[aln_len] = ref_aln[k];
                ref_read_aln_str->query_aln[aln_len] = read_aln[k];
                aln_len++;
            }
            i+=ref_del_len; j+=read_del_len;
            free(ref_aln);
        } else if (ref_cons_aln_str->query_aln[i] != 5 && cons_read_aln_str->target_aln[j] != 5) { // both are not gaps: match or mismatch
            ref_read_aln_str->target_aln[aln_len] = ref_cons_aln_str->target_aln[i];
            ref_read_aln_str->query_aln[aln_len] = cons_read_aln_str->query_aln[j];
            aln_len++; i++; j++;
        } else if (ref_cons_aln_str->query_aln[i] == 5 && cons_read_aln_str->target_aln[j] != 5) { // ref-cons is gap, cons-read is not gap
            ref_read_aln_str->target_aln[aln_len] = ref_cons_aln_str->target_aln[i];
            ref_read_aln_str->query_aln[aln_len] = 5;
            aln_len++; i++;
        } else if (ref_cons_aln_str->query_aln[i] != 5 && cons_read_aln_str->target_aln[j] == 5) { // ref-cons is not gap, cons-read is gap
            ref_read_aln_str->target_aln[aln_len] = 5;
            ref_read_aln_str->query_aln[aln_len] = cons_read_aln_str->query_aln[j];
            aln_len++; j++;
        } else {
            fprintf(stderr, "Error: ref_cons_aln_str->query_aln[%d]: %d cons_read_aln_str->target_aln[%d]: %d\n", i, ref_cons_aln_str->query_aln[i], j, cons_read_aln_str->target_aln[j]);
            exit(1);
        }
    }
    while (i < ref_cons_aln_str->aln_len) {
        ref_read_aln_str->target_aln[aln_len] = ref_cons_aln_str->target_aln[i];
        ref_read_aln_str->query_aln[aln_len] = 5;
        aln_len++; i++;
    }
    while (j < cons_read_aln_str->aln_len) {
        ref_read_aln_str->target_aln[aln_len] = 5;
        ref_read_aln_str->query_aln[aln_len] = cons_read_aln_str->query_aln[j];
        aln_len++; j++;
    }
    ref_read_aln_str->aln_len = aln_len;
    ref_read_aln_str->target_beg = ref_read_aln_str->target_end = ref_read_aln_str->query_beg = ref_read_aln_str->query_end = -1; // unset
    if (LONGCALLD_VERBOSE >= 3) { // print ref-read aln-str
        fprintf(stderr, "Ref-Cons AlnStr: len: %d\n", ref_cons_aln_str->aln_len);
        for (int k = 0; k < ref_cons_aln_str->aln_len; ++k) {
            fprintf(stderr, "%c", "ACGTN-"[ref_cons_aln_str->target_aln[k]]);
        } fprintf(stderr, "\n");
        for (int k = 0; k < ref_cons_aln_str->aln_len; ++k) {
            fprintf(stderr, "%c", "ACGTN-"[ref_cons_aln_str->query_aln[k]]);
        } fprintf(stderr, "\n");

        fprintf(stderr, "Cons-Read AlnStr: len: %d\n", cons_read_aln_str->aln_len);
        for (int k = 0; k < cons_read_aln_str->aln_len; ++k) {
            fprintf(stderr, "%c", "ACGTN-"[cons_read_aln_str->target_aln[k]]);
        } fprintf(stderr, "\n");
        for (int k = 0; k < cons_read_aln_str->aln_len; ++k) {
            fprintf(stderr, "%c", "ACGTN-"[cons_read_aln_str->query_aln[k]]);
        } fprintf(stderr, "\n");

        fprintf(stderr, "Ref-Read AlnStr: len: %d\n", aln_len);
        for (int k = 0; k < aln_len; ++k) {
            fprintf(stderr, "%c", "ACGTN-"[ref_read_aln_str->target_aln[k]]);
        } fprintf(stderr, "\n");
        for (int k = 0; k < aln_len; ++k) {
            fprintf(stderr, "%c", "ACGTN-"[ref_read_aln_str->query_aln[k]]);
        } fprintf(stderr, "\n");
    }
    return aln_len;
}

int wfa_collect_noisy_aln_str_no_ps_hap(const call_var_opt_t *opt, int n_reads, int *read_ids, int *lens, uint8_t **seqs, char **qnames, int *fully_covers,
                                        uint8_t *ref_seq, int ref_seq_len, int *clu_n_seqs, int **clu_read_ids, aln_str_t **aln_strs, int collect_ref_read_aln_str) {
    int *full_read_ids = (int*)malloc((n_reads+2) * sizeof(int));
    int *full_read_lens = (int*)malloc((n_reads+2) * sizeof(int));
    uint8_t **full_read_seqs = (uint8_t**)malloc((n_reads+2) * sizeof(uint8_t*));
    char **full_read_names = (char**)malloc((n_reads+2) * sizeof(char*));
    int *full_fully_covers = (int*)malloc((n_reads+2) * sizeof(int));
    int n_full_reads = 0;
    for (int i = 0; i < n_reads; ++i) {
        if (lens[i] <= 0 || LONGCALLD_NOISY_IS_BOTH_COVER(fully_covers[i]) == 0) continue;
        full_fully_covers[n_full_reads] = fully_covers[i];
        full_read_ids[n_full_reads] = i;
        full_read_lens[n_full_reads] = lens[i];
        full_read_seqs[n_full_reads] = seqs[i];
        full_read_names[n_full_reads] = qnames[i];
        n_full_reads++;
    }
    // abpoa
    int *cons_lens = (int*)malloc(2 * sizeof(int)); uint8_t **cons_seqs = (uint8_t**)malloc(2 * sizeof(uint8_t*));
    int *msa_seq_lens = (int*)malloc(2 * sizeof(int)); uint8_t ***msa_seqs = (uint8_t***)malloc(2 * sizeof(uint8_t**));
    for (int i = 0; i < 2; ++i) {
        cons_seqs[i] = NULL; cons_lens[i] = 0;
        msa_seqs[i] = (uint8_t**)malloc((n_reads+1) * sizeof(uint8_t*));
        for (int j = 0; j < n_reads+1; ++j) msa_seqs[i][j] = NULL;
        msa_seq_lens[i] = 0;
    }
    int n_cons = 0;
    if (n_full_reads == 0) goto collect_noisy_msa_cons_no_ps_hap_end;
    else {
        if (full_read_lens[0] >= opt->max_noisy_reg_len) goto collect_noisy_msa_cons_no_ps_hap_end;
    }

    n_cons = abpoa_aln_msa_cons(opt, n_full_reads, full_read_ids, full_read_seqs, full_read_lens, 2,
                                cons_lens, cons_seqs, clu_n_seqs, clu_read_ids, msa_seq_lens, msa_seqs);

    // re-do POA with ref_seq and cons
    for (int i = 0; i < n_cons; ++i) {
        aln_str_t *clu_aln_str = aln_strs[i];
        wfa_collect_aln_str(opt, ref_seq, ref_seq_len, cons_seqs[i], cons_lens[i], LONGCALLD_NOISY_BOTH_COVER, LONGCALLD_WFA_NO_HEURISTIC, LONGCALLD_WFA_AFFINE_2P, LONGCALLD_REF_CONS_ALN_STR(clu_aln_str));
        n_full_reads = 0;
        for (int j = 0; j < clu_n_seqs[i]; ++j) {
            int read_i = clu_read_ids[i][j];
            int read_id = read_ids[read_i];
            clu_read_ids[i][j] = read_id;
            // cons vs read
            // wfa_collect_aln_str(opt, cons_seqs[i], cons_lens[i], seqs[read_i], lens[read_i], fully_covers[read_i], LONGCALLD_WFA_NO_HEURISTIC, LONGCALLD_WFA_AFFINE_2P, LONGCALLD_CONS_READ_ALN_STR(clu_aln_str, n_full_reads));
            make_cons_read_aln_str(opt, msa_seqs[i][clu_n_seqs[i]], msa_seqs[i][j], msa_seq_lens[i], fully_covers[i], LONGCALLD_CONS_READ_ALN_STR(clu_aln_str, n_full_reads));
            if (collect_ref_read_aln_str)
                make_ref_read_aln_str(opt, LONGCALLD_REF_CONS_ALN_STR(clu_aln_str), LONGCALLD_CONS_READ_ALN_STR(clu_aln_str, n_full_reads), LONGCALLD_REF_READ_ALN_STR(clu_aln_str, n_full_reads));
            n_full_reads++;
        }
    } 
collect_noisy_msa_cons_no_ps_hap_end:
    free(full_read_lens); free(full_read_seqs); free(full_read_names); free(full_fully_covers); free(full_read_ids);
    for (int i = 0; i < 2; ++i) {
        for (int j = 0; j < n_reads+1; ++j) {
            if (msa_seqs[i][j] != NULL) free(msa_seqs[i][j]);
        } free(msa_seqs[i]);
    } free(msa_seqs); free(msa_seq_lens);
    for (int i = 0; i < 2; ++i) {
        if (cons_seqs[i] != NULL) free(cons_seqs[i]);
    } free(cons_lens); free(cons_seqs);
    return n_cons;
}

int ps_has_enough_reads(int ps_i, int **phase_set_to_hap_full_read_count, int **phase_set_to_hap_left_read_count, int **phase_set_to_hap_right_read_count,
                        int min_hap_full_read_count, int min_hap_all_read_count, int min_hap_partial_read_count) {
    if (ps_i < 0) return 0; // invalid phase set
    for (int hap=0; hap < 2; ++hap) {
        if ((phase_set_to_hap_full_read_count[ps_i][hap] < min_hap_full_read_count || phase_set_to_hap_left_read_count[ps_i][hap] + phase_set_to_hap_right_read_count[ps_i][hap] < min_hap_all_read_count)
          && (phase_set_to_hap_left_read_count[ps_i][hap] < min_hap_partial_read_count || phase_set_to_hap_right_read_count[ps_i][hap] < min_hap_partial_read_count))
          return 0;
    }
    return 1;
}

// collect phase set with both haps having >= min_minor_hap_read_count reads
hts_pos_t collect_phase_set_with_both_haps(int n_reads, int *read_haps, int *read_lens, hts_pos_t *phase_sets, int *fully_covers, int min_minor_hap_full_read_count, int min_minor_hap_all_read_count) {
    int n_uniq_phase_sets = 0, phase_set_i = 0;
    hts_pos_t *uniq_phase_sets = (hts_pos_t*)calloc(n_reads, sizeof(hts_pos_t));
    int **phase_set_to_hap_full_read_count = (int**)malloc(n_reads * sizeof(int*));
    int **phase_set_to_hap_all_read_count = (int**)malloc(n_reads * sizeof(int*)); // including clipped reads
    int **phase_set_to_hap_min_full_read_len = (int**)malloc(n_reads * sizeof(int*));
    for (int i = 0; i < n_reads; ++i) {
        phase_set_to_hap_full_read_count[i] = (int*)calloc(2, sizeof(int));
        phase_set_to_hap_all_read_count[i] = (int*)calloc(2, sizeof(int));
        phase_set_to_hap_min_full_read_len[i] = (int*)malloc(2 * sizeof(int));
        for (int j = 0; j < 2; ++j) phase_set_to_hap_min_full_read_len[i][j] = INT32_MAX;
    }
    // use one phase set if multiple phase sets exist, reads in other phase sets are considered as non-HAP
    for (int i = 0; i < n_reads; ++i) {
        if (read_haps[i] == 0) continue;
        phase_set_i = add_phase_set(phase_sets[i], uniq_phase_sets, &n_uniq_phase_sets);
        if (LONGCALLD_NOISY_IS_BOTH_COVER(fully_covers[i])) {
            phase_set_to_hap_full_read_count[phase_set_i][read_haps[i]-1]++;
            phase_set_to_hap_all_read_count[phase_set_i][read_haps[i]-1]++;
            if (phase_set_to_hap_min_full_read_len[phase_set_i][read_haps[i]-1] > read_lens[i]) 
                phase_set_to_hap_min_full_read_len[phase_set_i][read_haps[i]-1] = read_lens[i];
        } else if (LONGCALLD_NOISY_IS_LEFT_COVER(fully_covers[i]) || LONGCALLD_NOISY_IS_RIGHT_COVER(fully_covers[i])) {
            if (read_lens[i] >= phase_set_to_hap_min_full_read_len[phase_set_i][read_haps[i]-1]) {
                phase_set_to_hap_all_read_count[phase_set_i][read_haps[i]-1]++;
            }
        }
    }
    hts_pos_t max_ps = -1; int max_ps_i = -1; 
    int max_ps_full_read_count1 = -1, max_ps_full_read_count2 = -1;
    for (int i = 0; i < n_uniq_phase_sets; ++i) {
        int phase_set_full_read_count1 = phase_set_to_hap_full_read_count[i][0] < phase_set_to_hap_full_read_count[i][1] ? phase_set_to_hap_full_read_count[i][0] : phase_set_to_hap_full_read_count[i][1];
        int phase_set_full_read_count2 = phase_set_to_hap_full_read_count[i][0] > phase_set_to_hap_full_read_count[i][1] ? phase_set_to_hap_full_read_count[i][0] : phase_set_to_hap_full_read_count[i][1];
        if (phase_set_full_read_count1 > max_ps_full_read_count1) {
            max_ps_full_read_count1 = phase_set_full_read_count1;
            max_ps_full_read_count2 = phase_set_full_read_count2;
            max_ps = uniq_phase_sets[i];
            max_ps_i = i;
        } else if (phase_set_full_read_count1 == max_ps_full_read_count1 && phase_set_full_read_count2 > max_ps_full_read_count2) {
            max_ps_full_read_count2 = phase_set_full_read_count2;
            max_ps = uniq_phase_sets[i];
            max_ps_i = i;
        }
    }
    if (max_ps_full_read_count1 < min_minor_hap_full_read_count) max_ps = -1;
    if (max_ps != -1 && max_ps_i != -1) {
        if (phase_set_to_hap_all_read_count[max_ps_i][0] < min_minor_hap_all_read_count ||
            phase_set_to_hap_all_read_count[max_ps_i][1] < min_minor_hap_all_read_count) max_ps = -1;
    }
    for (int i = 0; i < n_reads; ++i) {
        free(phase_set_to_hap_full_read_count[i]);
        free(phase_set_to_hap_all_read_count[i]);
        free(phase_set_to_hap_min_full_read_len[i]);
    } free(uniq_phase_sets); free(phase_set_to_hap_full_read_count); free(phase_set_to_hap_all_read_count); free(phase_set_to_hap_min_full_read_len);
    return max_ps;
}

// pairwise alignment of cons and ref to append ref to the MSA
// total: 1+ n_reads*2
// 1. ref vs cons: 1
// 2. cons vs n_reads: n_reads
// 3. ref vs n_reads: n_reads
int wfa_collect_noisy_aln_str_with_ps_hap(const call_var_opt_t *opt, int sampling_reads, int n_reads, int *noisy_read_ids, int *lens, uint8_t **seqs, uint8_t *strands, uint8_t **quals, char **names,
                                          int *haps, hts_pos_t *phase_sets, int *fully_covers, hts_pos_t ps, int min_hap_full_reads, int min_hap_all_reads, uint8_t *ref_seq, int ref_seq_len,
                                          int *clu_n_seqs, int **clu_read_ids, aln_str_t **aln_strs, int collect_ref_read_aln_str) {
    // given specific phase_set, collect consensus sequences for each haplotype
    // optional: for long noisy regions, skip read with large edit distance to first read
    int n_cons = 0;
    int total_n_reads = n_reads+2, n_ps_hap_reads = 0;
    int *ps_hap_read_ids = (int*)malloc(total_n_reads * sizeof(int));
    int *ps_hap_read_lens = (int*)malloc(total_n_reads * sizeof(int));
    uint8_t **ps_hap_read_seqs = (uint8_t**)malloc(total_n_reads * sizeof(uint8_t*));
    uint8_t *ps_hap_read_strands = (uint8_t*)malloc(total_n_reads * sizeof(uint8_t));
    uint8_t **ps_hap_read_quals = (uint8_t**)malloc(total_n_reads * sizeof(uint8_t*));
    char **ps_hap_read_names = (char**)malloc(total_n_reads * sizeof(char*));
    int *ps_hap_full_covers = (int*)calloc(total_n_reads, sizeof(int));
    int *cons_lens = (int*)malloc(2 * sizeof(int)); uint8_t **cons_seqs = (uint8_t**)malloc(2 * sizeof(uint8_t*));
    int *msa_seq_lens = (int*)malloc(2 * sizeof(int)); uint8_t ***msa_seqs = (uint8_t***)malloc(2 * sizeof(uint8_t**));
    for (int i = 0; i < 2; ++i) {
        cons_seqs[i] = NULL; cons_lens[i] = 0;
        msa_seqs[i] = (uint8_t**)malloc((n_reads+1) * sizeof(uint8_t*));
        for (int j = 0; j < n_reads+1; ++j) msa_seqs[i][j] = NULL;
        msa_seq_lens[i] = 0;
    }

    int hp_flank_start, hp_flank_end, hp_len, use_non_full = 1;
    if (is_homopolymer(ref_seq, ref_seq_len, opt->noisy_reg_flank_len, &hp_flank_start, &hp_flank_end, &hp_len)) {
        use_non_full = 0;
    }
    for (int hap=1; hap<=2; ++hap) {
        // check if we have enough full-cover reads
        n_ps_hap_reads = 0;
        for (int i = 0; i < n_reads; ++i) {
            if (lens[i] <= 0 || phase_sets[i] != ps || haps[i] != hap) continue;
            if (use_non_full == 0 && LONGCALLD_NOISY_IS_BOTH_COVER(fully_covers[i]) == 0) continue;
            ps_hap_read_ids[n_ps_hap_reads] = noisy_read_ids[i];
            ps_hap_read_lens[n_ps_hap_reads] = lens[i];
            ps_hap_read_seqs[n_ps_hap_reads] = seqs[i];
            ps_hap_read_strands[n_ps_hap_reads] = strands[i];
            ps_hap_read_quals[n_ps_hap_reads] = quals[i];
            ps_hap_full_covers[n_ps_hap_reads] = fully_covers[i];
            ps_hap_read_names[n_ps_hap_reads] = names[i];
            n_ps_hap_reads++;
        }
        if (ps_hap_read_lens[0] >= opt->max_noisy_reg_len) {
            if (LONGCALLD_VERBOSE >= 2) fprintf(stderr, "SkipRegion: %" PRIi64 " %d %d %d\n", ps, hap, n_ps_hap_reads, ps_hap_read_lens[0]);
            break;
        }
        if (LONGCALLD_VERBOSE >= 2) fprintf(stderr, "PS: %" PRIi64 " HAP: %d n_reads: %d\n", ps, hap, n_ps_hap_reads);
        if (n_ps_hap_reads == 0) continue;
        // collect consensus sequences
        n_cons += abpoa_partial_aln_msa_cons(opt, NULL, sampling_reads, n_ps_hap_reads, ps_hap_read_ids, ps_hap_read_seqs, ps_hap_read_quals, ps_hap_read_lens, ps_hap_full_covers, ps_hap_read_names,
                                             1, cons_lens+hap-1, cons_seqs+hap-1, clu_n_seqs+hap-1, clu_read_ids+hap-1, msa_seq_lens+hap-1, msa_seqs[hap-1]);
    }

    if (n_cons != 2) n_cons = 0;
    else { // collect aln_strs
        for (int hap=1; hap<=2; ++hap) {
            // ref vs cons
            aln_str_t *clu_aln_str = aln_strs[hap-1];
            // fprintf(stderr, "WFA cons-ref align for HAP: %d %d vs %d\n", hap, ref_seq_len, cons_lens[hap-1]);
            wfa_collect_aln_str(opt, ref_seq, ref_seq_len, cons_seqs[hap-1], cons_lens[hap-1], LONGCALLD_NOISY_BOTH_COVER, LONGCALLD_WFA_NO_HEURISTIC, LONGCALLD_WFA_AFFINE_2P, LONGCALLD_REF_CONS_ALN_STR(clu_aln_str));
            n_ps_hap_reads = 0;
            for (int i = 0; i < n_reads; ++i) {
                if (lens[i] <= 0 || phase_sets[i] != ps || haps[i] != hap) continue;
                if (use_non_full == 0 && LONGCALLD_NOISY_IS_BOTH_COVER(fully_covers[i]) == 0) continue;
                // cons vs read XXX make?
                // heuristic = 1; affine_2p = 0; // 
                // wfa_collect_aln_str(opt, cons_seqs[hap-1], cons_lens[hap-1], seqs[i], lens[i], fully_covers[i], LONGCALLD_WFA_NO_HEURISTIC, LONGCALLD_WFA_AFFINE_2P, LONGCALLD_CONS_READ_ALN_STR(clu_aln_str, n_ps_hap_reads));
                make_cons_read_aln_str(opt, msa_seqs[hap-1][clu_n_seqs[hap-1]], msa_seqs[hap-1][n_ps_hap_reads], msa_seq_lens[hap-1], fully_covers[i], LONGCALLD_CONS_READ_ALN_STR(clu_aln_str, n_ps_hap_reads));
                // ref vs read
                if (collect_ref_read_aln_str)
                    // fprintf(stderr, "Make ref-read aln str for HAP: %d read_id: %d, %d vs %d (%d)\n", hap, noisy_read_ids[i], lens[i], ref_seq_len, fully_covers[i]);
                    make_ref_read_aln_str(opt, LONGCALLD_REF_CONS_ALN_STR(clu_aln_str), LONGCALLD_CONS_READ_ALN_STR(clu_aln_str, n_ps_hap_reads), LONGCALLD_REF_READ_ALN_STR(clu_aln_str, n_ps_hap_reads));
                n_ps_hap_reads++;
            }
            if (LONGCALLD_VERBOSE >=2 ) fprintf(stderr, "With Ref+Cons PS: %" PRIi64 " HAP: %d n_reads: %d\n", ps, hap, n_ps_hap_reads);
            if (n_ps_hap_reads == 0) continue;
        }
    }
    free(ps_hap_read_ids); free(ps_hap_read_lens); free(ps_hap_read_seqs); free(ps_hap_read_strands); free(ps_hap_read_quals); free(ps_hap_full_covers); free(ps_hap_read_names);
    for (int i = 0; i < 2; ++i) {
        if (cons_seqs[i] != NULL) free(cons_seqs[i]);
        if (msa_seqs[i] != NULL) {
            for (int j = 0; j < n_reads+1; ++j) {
                if (msa_seqs[i][j] != NULL) free(msa_seqs[i][j]);
            }
            free(msa_seqs[i]);
        }
    } free(cons_lens); free(cons_seqs); free(msa_seq_lens); free(msa_seqs);
    return n_cons;
}

int collect_noisy_read_info(const call_var_opt_t *opt, bam_chunk_t *chunk, hts_pos_t reg_beg, hts_pos_t reg_end, int noisy_reg_i, int n_noisy_reg_reads, int *noisy_reg_reads, int **read_lens,
                            uint8_t ***read_seqs, uint8_t **strands, uint8_t ***read_quals, char ***read_names, int **fully_covers, int **read_id_to_full_covers, 
                            int **read_reg_beg, int **read_reg_end, int **read_haps, hts_pos_t **phase_sets) {
    *read_lens = (int*)calloc(n_noisy_reg_reads, sizeof(int));
    *read_seqs = (uint8_t**)malloc(n_noisy_reg_reads * sizeof(uint8_t*));
    *strands = (uint8_t*)malloc(n_noisy_reg_reads * sizeof(uint8_t));
    *read_quals = (uint8_t**)malloc(n_noisy_reg_reads * sizeof(uint8_t*));
    *read_names = (char**)malloc(n_noisy_reg_reads * sizeof(char*));
    *fully_covers = (int*)calloc(n_noisy_reg_reads, sizeof(int));
    *read_id_to_full_covers = (int*)calloc(chunk->n_reads, sizeof(int));
    *read_reg_beg = (int*)calloc(chunk->n_reads, sizeof(int));
    *read_reg_end = (int*)calloc(chunk->n_reads, sizeof(int));
    *read_haps = (int*)calloc(n_noisy_reg_reads, sizeof(int));
    *phase_sets = (hts_pos_t*)calloc(n_noisy_reg_reads, sizeof(hts_pos_t));

    for (int i = 0; i < n_noisy_reg_reads; ++i) {
        int read_id = noisy_reg_reads[i];
        digar_t *read_digars = chunk->digars+read_id; int n_digar = read_digars->n_digar; digar1_t *digars = read_digars->digars;
        hts_pos_t reg_digar_beg = -1, reg_digar_end = -1;
        int reg_read_beg = 0, reg_read_end = digar2qlen(read_digars)-1;
        if (read_digars->digars[0].type == BAM_CHARD_CLIP) reg_read_beg = read_digars->digars[0].len;
        if (read_digars->digars[n_digar-1].type == BAM_CHARD_CLIP) reg_read_end = read_digars->digars[n_digar-1].qi - 1;
        if (LONGCALLD_VERBOSE >= 2) (*read_names)[i] = bam_get_qname(chunk->reads[read_id]);
        else (*read_names)[i] = NULL;
        (*strands)[i] = read_digars->is_rev;
        int beg_is_del = 0, end_is_del = 0, cover = 0;
        for (int digar_i = 0; digar_i < n_digar; ++digar_i) {
            hts_pos_t digar_beg = digars[digar_i].pos, digar_end;
            int op = digars[digar_i].type, len = digars[digar_i].len, qi = digars[digar_i].qi;
            if (op == BAM_CSOFT_CLIP || op == BAM_CHARD_CLIP) continue;
            if (op == BAM_CDIFF || op == BAM_CEQUAL || op == BAM_CDEL) digar_end = digar_beg + len - 1;
            else digar_end = digar_beg;
            if (digar_beg > reg_end) break;
            if (digar_end < reg_beg) continue;
            if (digar_beg <= reg_beg && digar_end >= reg_beg) {
                if (op == BAM_CDEL) {
                    reg_digar_beg = reg_beg;
                    reg_read_beg = qi; // qi is on the right side of the DEL
                    if (LONGCALLD_VERBOSE >= 2) fprintf(stderr, "%s noisyReg left boundary is DEL %s:%" PRId64 "-%" PRId64 "\n", (*read_names)[i], chunk->tname, reg_beg, reg_end);
                    if (len > opt->noisy_reg_flank_len) beg_is_del = 1;
                } else {
                    reg_digar_beg = reg_beg;
                    reg_read_beg = qi + (reg_beg - digar_beg);
                }
            }
            if (digar_beg <= reg_end && digar_end >= reg_end) {
                if (op == BAM_CDEL) {
                    reg_digar_end = reg_end;
                    reg_read_end = qi-1; // qi is the one the left side of the DEL
                    if (LONGCALLD_VERBOSE >= 2) fprintf(stderr, "%s noisyReg right boundary is DEL %s:%" PRId64 "-%" PRId64 "\n", (*read_names)[i], chunk->tname, reg_beg, reg_end);
                    if (len > opt->noisy_reg_flank_len) end_is_del = 1;
                } else {
                    reg_digar_end = reg_end;
                    reg_read_end = qi + (reg_end - digar_beg);
                }
            }
        }
        if (reg_digar_beg == reg_beg && reg_digar_end == reg_end) {
            if (beg_is_del == 0 && end_is_del == 0) cover = LONGCALLD_NOISY_LEFT_COVER | LONGCALLD_NOISY_RIGHT_COVER; // (*fully_covers)[i] = 3;
            else if (beg_is_del == 0 && end_is_del == 1) cover = LONGCALLD_NOISY_LEFT_COVER | LONGCALLD_NOISY_RIGHT_GAP; // (*fully_covers)[i] = 1;
            else if (beg_is_del == 1 && end_is_del == 0) cover = LONGCALLD_NOISY_LEFT_GAP | LONGCALLD_NOISY_RIGHT_COVER; // (*fully_covers)[i] = 2;
            else cover = LONGCALLD_NOISY_LEFT_GAP | LONGCALLD_NOISY_RIGHT_GAP; // (*fully_covers)[i] = 0;
        } else if (reg_digar_beg == reg_beg) {
            if (beg_is_del) cover = LONGCALLD_NOISY_LEFT_GAP; // (*fully_covers)[i] = 0;
            else cover = LONGCALLD_NOISY_LEFT_COVER; // (*fully_covers)[i] = 1;
        } else if (reg_digar_end == reg_end) {
            if (end_is_del) cover = LONGCALLD_NOISY_RIGHT_GAP; // (*fully_covers)[i] = 0;
            else cover = LONGCALLD_NOISY_RIGHT_COVER; // (*fully_covers)[i] = 2;
        } else cover = 0; // (*fully_covers)[i] = 0;
        // if (2*(reg_read_end-reg_read_beg+1) < (reg_end-reg_beg+1)) return 0;
        (*read_seqs)[i] = (uint8_t*)malloc((reg_read_end - reg_read_beg + 1) * sizeof(uint8_t));
        (*read_quals)[i] = (uint8_t*)malloc((reg_read_end - reg_read_beg + 1) * sizeof(uint8_t));
        for (int j = reg_read_beg; j <= reg_read_end; ++j) {
            (*read_seqs)[i][j-reg_read_beg] = seq_nt16_int[bam_seqi(read_digars->bseq, j)];
            (*read_quals)[i][j-reg_read_beg] = read_digars->qual[j];
        }
        (*read_lens)[i] = reg_read_end - reg_read_beg + 1;
        (*read_haps)[i] = chunk->haps[read_id];
        (*phase_sets)[i] = chunk->phase_sets[read_id];
        (*fully_covers)[i] = cover;
        (*read_id_to_full_covers)[read_id] = cover;
        (*read_reg_beg)[read_id] = reg_read_beg; (*read_reg_end)[read_id] = reg_read_end;
    }
    return 0;
}

digar1_t *collect_left_digars(digar_t *digar, int read_noisy_beg, hts_pos_t ref_noisy_beg, int *n_left_digars) {
    digar1_t *left_digars = NULL;
    *n_left_digars = 0; int m_left_digar = 0;
    for (int i = 0; i < digar->n_digar; ++i) {
        int op = digar->digars[i].type, qi = digar->digars[i].qi, digar_qi_end; 
        hts_pos_t digar_ref_beg = digar->digars[i].pos, digar_ref_end;
        // collect left clipping
        if (i == 0 && (op == BAM_CSOFT_CLIP || op == BAM_CHARD_CLIP)) {
            left_digars = push_digar0(left_digars, n_left_digars, &m_left_digar, digar->digars[i]);
            continue;
        }

        if (op == BAM_CDIFF || op == BAM_CEQUAL || op == BAM_CINS) digar_qi_end = qi + digar->digars[i].len - 1;
        else digar_qi_end = qi;
        if (op == BAM_CDIFF || op == BAM_CEQUAL || op == BAM_CDEL) digar_ref_end = digar->digars[i].pos + digar->digars[i].len - 1;
        else digar_ref_end = digar->digars[i].pos;
        if (qi >= read_noisy_beg && digar_ref_beg >= ref_noisy_beg) break; // reach noisy region

        if (digar_qi_end < read_noisy_beg && digar_ref_end < ref_noisy_beg) // full digar
            left_digars = push_digar0(left_digars, n_left_digars, &m_left_digar, digar->digars[i]);
        else if (digar_qi_end >= read_noisy_beg || digar_ref_end >= ref_noisy_beg) { // partial digar, chop X=ID
            if (op == BAM_CINS || op == BAM_CEQUAL || op == BAM_CDIFF) {
                digar1_t d = digar->digars[i]; d.len = read_noisy_beg - d.qi; // chop digar
                left_digars = push_digar0(left_digars, n_left_digars, &m_left_digar, d);
            } else if (op == BAM_CDEL) {
                digar1_t d = digar->digars[i]; d.len = ref_noisy_beg - d.pos; // chop digar
                left_digars = push_digar0(left_digars, n_left_digars, &m_left_digar, d);
            }
            break;
        } else {
            fprintf(stderr, "Error: digar_ref_end: %" PRId64 " ref_noisy_beg: %" PRId64 "\n", digar_ref_end, ref_noisy_beg);
            exit(1);
        }
    }
    return left_digars;
}

digar1_t *collect_right_digars(digar_t *digar, int read_noisy_end, hts_pos_t ref_noisy_end, int *n_right_digars) {
    digar1_t *right_digars = NULL; *n_right_digars = 0; int m_right_digar = 0;
    for (int i = 0; i < digar->n_digar; ++i) {
        int qi = digar->digars[i].qi, op = digar->digars[i].type, digar_qi_end;
        hts_pos_t digar_ref_beg = digar->digars[i].pos, digar_ref_end;
        // collect right clipping
        if (i == digar->n_digar-1 && (op == BAM_CSOFT_CLIP || op == BAM_CHARD_CLIP)) {
            right_digars = push_digar0(right_digars, n_right_digars, &m_right_digar, digar->digars[i]);
            continue;
        }
        if (op == BAM_CDIFF || op == BAM_CEQUAL || op == BAM_CINS) digar_qi_end = qi + digar->digars[i].len - 1;
        else digar_qi_end = qi;
        if (op == BAM_CDIFF || op == BAM_CEQUAL || op == BAM_CDEL) digar_ref_end = digar->digars[i].pos + digar->digars[i].len - 1;
        else digar_ref_end = digar->digars[i].pos;
        if (digar_qi_end <= read_noisy_end && digar_ref_end <= ref_noisy_end) continue; // left side of noisy region

        if (qi > read_noisy_end && digar_ref_beg > ref_noisy_end) // full digar
            right_digars = push_digar0(right_digars, n_right_digars, &m_right_digar, digar->digars[i]);
        else if (qi <= read_noisy_end || digar_ref_beg <= ref_noisy_end) { // partial digar, chop IDX
            if (op == BAM_CINS || op == BAM_CEQUAL || op == BAM_CDIFF) {
                digar1_t d = digar->digars[i]; d.len = digar_qi_end - read_noisy_end; // chop digar
                d.qi = read_noisy_end + 1; // chop qi
                if (op == BAM_CINS) {
                    uint8_t *new_alt_seq = (uint8_t*)malloc(d.len * sizeof(uint8_t));
                    for (int j = 0; j < d.len; ++j) {
                        new_alt_seq[j] = digar->digars[i].alt_seq[j+read_noisy_end+1-qi];
                    }
                    free(d.alt_seq); d.alt_seq = new_alt_seq; // chop alt_seq
                } else d.pos = ref_noisy_end + 1; // chop pos
                right_digars = push_digar0(right_digars, n_right_digars, &m_right_digar, d);
            } else if (op == BAM_CDEL) {
                digar1_t d = digar->digars[i]; d.len = digar_ref_end - ref_noisy_end; // chop digar
                d.pos = ref_noisy_end + 1; // chop pos
                right_digars = push_digar0(right_digars, n_right_digars, &m_right_digar, d);
            }
        } else {
            fprintf(stderr, "Error: digar_ref_end: %" PRId64 " ref_noisy_end: %" PRId64 "\n", digar_ref_end, ref_noisy_end);
            exit(1);
        }
    }
    return right_digars;
}

digar1_t *collect_full_msa_digars(int read_beg, int read_end, hts_pos_t ref_beg, int msa_len, uint8_t *read_str, uint8_t *ref_str, int *n_msa_digars) {
    digar1_t *msa_digars = NULL;
    *n_msa_digars = 0; int m_msa_digars = 0;
    if (msa_len <= 0) return NULL;
    int read_pos = read_beg; hts_pos_t ref_pos = ref_beg;
    int left_read_start = 0, right_read_end = msa_len - 1;
    for (int i = 0; i < msa_len; ++i) {
        if (read_str[i] == 5 && ref_str[i] == 5) continue; // both are non-ACGTN
        if (read_str[i] != 5 && ref_str[i] != 5) {
            if (i >= left_read_start && i <= right_read_end) {
                if (read_str[i] == ref_str[i]) { // =
                    digar1_t d = (digar1_t){.pos = ref_pos, .qi = read_pos, .type = BAM_CEQUAL, .len = 1, .alt_seq = NULL, .is_low_qual = 0};
                    msa_digars = push_digar0(msa_digars, n_msa_digars, &m_msa_digars, d);
                } else { // X
                    digar1_t d = (digar1_t){.pos = ref_pos, .qi = read_pos, .type = BAM_CDIFF, .len = 1, .alt_seq = NULL, .is_low_qual = 0};
                    d.alt_seq = (uint8_t*)malloc(1 * sizeof(uint8_t));
                    d.alt_seq[0] = read_str[i];
                    msa_digars = push_digar0(msa_digars, n_msa_digars, &m_msa_digars, d);
                }
            }
            read_pos++; ref_pos++;
        } else if (read_str[i] != 5) { // INS: ref_str[i] == 5
            if (i >= left_read_start && i <= right_read_end) {
                digar1_t d = (digar1_t){.pos = ref_pos, .qi = read_pos, .type = BAM_CINS, .len = 1, .alt_seq = NULL, .is_low_qual = 0};
                d.alt_seq = (uint8_t*)malloc(1 * sizeof(uint8_t));
                d.alt_seq[0] = read_str[i];
                msa_digars = push_digar0(msa_digars, n_msa_digars, &m_msa_digars, d);
            }
            read_pos++;
        } else { // DEL read_str[i] == 5, ref_str[i] != 5
            if (i >= left_read_start && i <= right_read_end) {
                digar1_t d = (digar1_t){.pos = ref_pos, .qi = read_pos, .type = BAM_CDEL, .len = 1, .alt_seq = NULL, .is_low_qual = 0};
                msa_digars = push_digar0(msa_digars, n_msa_digars, &m_msa_digars, d);
            }
            ref_pos++;
        }
    }
    return msa_digars;
}

digar1_t *collect_left_msa_digars(int read_beg, int read_end, int qlen, hts_pos_t ref_beg, int msa_len, uint8_t *read_str, uint8_t *ref_str, int *n_msa_digars) {
    digar1_t *msa_digars = NULL;
    *n_msa_digars = 0; int m_msa_digars = 0;
    if (msa_len <= 0) return NULL;
    int read_pos = read_beg, read_end_pos = read_end; hts_pos_t ref_pos = ref_beg;
    int left_read_start = 0, right_read_end = msa_len - 1, right_skipped_read_base = 0;
    // set read_end_pos right_read_end
    read_end_pos = read_pos-1; int is_covered_by_ref = 0;
    for (int i = msa_len-1; i >= 0; --i) {
        if (ref_str[i] != 5) is_covered_by_ref = 1;
        if (is_covered_by_ref && read_str[i] != 5) {
            right_read_end = i; // find the last non-5 base
            break;
        } else if (!is_covered_by_ref && read_str[i] != 5) {
            right_skipped_read_base++;
        }
    }
    for (int i = 0; i < msa_len; ++i) {
        if (read_str[i] != 5) read_end_pos++;
    }
    for (int i = 0; i < msa_len; ++i) {
        if (read_str[i] == 5 && ref_str[i] == 5) continue; // both are non-ACGTN
        if (read_str[i] != 5 && ref_str[i] != 5) {
            if (i >= left_read_start && i <= right_read_end) {
                if (read_str[i] == ref_str[i]) { // =
                    digar1_t d = (digar1_t){.pos = ref_pos, .qi = read_pos, .type = BAM_CEQUAL, .len = 1, .alt_seq = NULL, .is_low_qual = 0};
                    msa_digars = push_digar0(msa_digars, n_msa_digars, &m_msa_digars, d);
                } else { // X
                    digar1_t d = (digar1_t){.pos = ref_pos, .qi = read_pos, .type = BAM_CDIFF, .len = 1, .alt_seq = NULL, .is_low_qual = 0};
                    d.alt_seq = (uint8_t*)malloc(1 * sizeof(uint8_t));
                    d.alt_seq[0] = read_str[i];
                    msa_digars = push_digar0(msa_digars, n_msa_digars, &m_msa_digars, d);
                }
            }
            read_pos++; ref_pos++;
        } else if (read_str[i] != 5) { // INS: ref_str[i] == 5
            if (i >= left_read_start && i <= right_read_end) {
                digar1_t d = (digar1_t){.pos = ref_pos, .qi = read_pos, .type = BAM_CINS, .len = 1, .alt_seq = NULL, .is_low_qual = 0};
                d.alt_seq = (uint8_t*)malloc(1 * sizeof(uint8_t));
                d.alt_seq[0] = read_str[i];
                msa_digars = push_digar0(msa_digars, n_msa_digars, &m_msa_digars, d);
            }
            read_pos++;
        } else { // DEL read_str[i] == 5, ref_str[i] != 5
            if (i >= left_read_start && i <= right_read_end) {
                digar1_t d = (digar1_t){.pos = ref_pos, .qi = read_pos, .type = BAM_CDEL, .len = 1, .alt_seq = NULL, .is_low_qual = 0};
                msa_digars = push_digar0(msa_digars, n_msa_digars, &m_msa_digars, d);
            }
            ref_pos++;
        }
    }
    if (read_end_pos < qlen-1 || right_skipped_read_base > 0) {
        // push back S digar XXX do extend alignment
        digar1_t d = (digar1_t){.pos = ref_pos, .qi = read_end_pos+1, .type = BAM_CSOFT_CLIP, .len = qlen-1-read_end_pos + right_skipped_read_base, .alt_seq = NULL, .is_low_qual = 0};
        msa_digars = push_digar0(msa_digars, n_msa_digars, &m_msa_digars, d);
    }
    return msa_digars;
}

digar1_t *collect_right_msa_digars(int read_beg, int read_end, hts_pos_t ref_beg, hts_pos_t ref_end, int msa_len, uint8_t *read_str, uint8_t *ref_str, int *n_msa_digars) {
    digar1_t *msa_digars = NULL;
    *n_msa_digars = 0; int m_msa_digars = 0;
    if (msa_len <= 0) return NULL;
    int read_pos = read_beg; hts_pos_t _ref_pos, ref_pos;
    int left_read_start = 0, right_read_end = msa_len - 1, left_skipped_read_base = 0;
        // set read_pos && left_read_start
    read_pos = read_end+1; _ref_pos = ref_end+1; ref_pos = ref_beg;
    int is_covered_by_ref = 0;
    for (int i = 0; i < msa_len; ++i) {
        if (ref_str[i] != 5) is_covered_by_ref = 1;
        if (is_covered_by_ref && read_str[i] != 5) {
            left_read_start = i; // find the first non-5 base
            break;
        } else if (!is_covered_by_ref && read_str[i] != 5) {
            left_skipped_read_base++;
        }
    }
    for (int i = msa_len-1; i>=0; --i) {
        if (ref_str[i] != 5) _ref_pos--;
        if (read_str[i] != 5) {
            read_pos--;
            ref_pos = _ref_pos;
        }
    }
    if (read_pos > 0 || left_skipped_read_base > 0) { // push front S digar XXX do extend alignment
        digar1_t d = (digar1_t){.pos = ref_pos, .qi = 0, .type = BAM_CSOFT_CLIP, .len = read_pos+left_skipped_read_base, .alt_seq = NULL, .is_low_qual = 0};
        msa_digars = push_digar0(msa_digars, n_msa_digars, &m_msa_digars, d);
    }
    read_pos += left_skipped_read_base;
    // for (int i = 0; i < msa_len; ++i) {
    for (int i = left_read_start; i <= right_read_end; ++i) {
        if (read_str[i] == 5 && ref_str[i] == 5) continue; // both are non-ACGTN
        if (read_str[i] != 5 && ref_str[i] != 5) {
            if (read_str[i] == ref_str[i]) { // =
                digar1_t d = (digar1_t){.pos = ref_pos, .qi = read_pos, .type = BAM_CEQUAL, .len = 1, .alt_seq = NULL, .is_low_qual = 0};
                msa_digars = push_digar0(msa_digars, n_msa_digars, &m_msa_digars, d);
            } else { // X
                digar1_t d = (digar1_t){.pos = ref_pos, .qi = read_pos, .type = BAM_CDIFF, .len = 1, .alt_seq = NULL, .is_low_qual = 0};
                d.alt_seq = (uint8_t*)malloc(1 * sizeof(uint8_t));
                d.alt_seq[0] = read_str[i];
                msa_digars = push_digar0(msa_digars, n_msa_digars, &m_msa_digars, d);
            }
            read_pos++; ref_pos++;
        } else if (read_str[i] != 5) { // INS: ref_str[i] == 5
            digar1_t d = (digar1_t){.pos = ref_pos, .qi = read_pos, .type = BAM_CINS, .len = 1, .alt_seq = NULL, .is_low_qual = 0};
            d.alt_seq = (uint8_t*)malloc(1 * sizeof(uint8_t));
            d.alt_seq[0] = read_str[i];
            msa_digars = push_digar0(msa_digars, n_msa_digars, &m_msa_digars, d);
            read_pos++;
        } else { // DEL read_str[i] == 5, ref_str[i] != 5
            digar1_t d = (digar1_t){.pos = ref_pos, .qi = read_pos, .type = BAM_CDEL, .len = 1, .alt_seq = NULL, .is_low_qual = 0};
            msa_digars = push_digar0(msa_digars, n_msa_digars, &m_msa_digars, d);
            ref_pos++;
        }
    }
    return msa_digars;
}

void update_digars_from_msa1(char *tname, digar_t *digar, int msa_len, uint8_t *ref_str, uint8_t *read_str, int full_cover, hts_pos_t noisy_reg_beg, hts_pos_t noisy_reg_end, int read_beg, int read_end) {
    int new_n_digars = 0, new_m_digar = 0; digar1_t *new_digars = NULL;
    int old_left_n_digars, old_right_n_digars, msa_n_digars;
    digar1_t *old_left_digars=NULL, *old_right_digars=NULL, *msa_digars=NULL;
    if (LONGCALLD_NOISY_IS_NOT_COVER(full_cover)) return;
    else if (LONGCALLD_NOISY_IS_BOTH_COVER(full_cover) ||  // both left & right
             (LONGCALLD_NOISY_IS_LEFT_COVER(full_cover) && LONGCALLD_NOISY_IS_RIGHT_GAP(full_cover)) ||
             (LONGCALLD_NOISY_IS_RIGHT_COVER(full_cover) && LONGCALLD_NOISY_IS_LEFT_GAP(full_cover))) {
        // chop digar
        old_left_digars = collect_left_digars(digar, read_beg, noisy_reg_beg, &old_left_n_digars);
        old_right_digars = collect_right_digars(digar, read_end, noisy_reg_end, &old_right_n_digars);
        msa_digars = collect_full_msa_digars(read_beg, read_end, noisy_reg_beg, msa_len, read_str, ref_str, &msa_n_digars);
        // push to new digar
        for (int i = 0; i < old_left_n_digars; ++i) new_digars = push_digar_alt_seq(new_digars, &new_n_digars, &new_m_digar, old_left_digars[i]);
        for (int i = 0; i < msa_n_digars; ++i) new_digars = push_digar0(new_digars, &new_n_digars, &new_m_digar, msa_digars[i]);
        for (int i = 0; i < old_right_n_digars; ++i) new_digars = push_digar_alt_seq(new_digars, &new_n_digars, &new_m_digar, old_right_digars[i]);
    } else if (LONGCALLD_NOISY_IS_LEFT_COVER(full_cover)) { // left-cover: update digar
        // chop digar
        old_left_digars = collect_left_digars(digar, read_beg, noisy_reg_beg, &old_left_n_digars);
        msa_digars = collect_left_msa_digars(read_beg, read_end, digar->qlen, noisy_reg_beg, msa_len, read_str, ref_str, &msa_n_digars);
        // push to new digar
        for (int i = 0; i < old_left_n_digars; ++i) new_digars = push_digar_alt_seq(new_digars, &new_n_digars, &new_m_digar, old_left_digars[i]);
        for (int i = 0; i < msa_n_digars; ++i) new_digars = push_digar0(new_digars, &new_n_digars, &new_m_digar, msa_digars[i]);
    } else if (LONGCALLD_NOISY_IS_RIGHT_COVER(full_cover)) { // right-cover: update digar and start pos
        // chop digar
        old_right_digars = collect_right_digars(digar, read_end, noisy_reg_end, &old_right_n_digars);
        msa_digars = collect_right_msa_digars(read_beg, read_end, noisy_reg_beg, noisy_reg_end, msa_len, read_str, ref_str, &msa_n_digars);
        // push to new digar
        for (int i = 0; i < msa_n_digars; ++i) new_digars = push_digar0(new_digars, &new_n_digars, &new_m_digar, msa_digars[i]);
        for (int i = 0; i < old_right_n_digars; ++i) new_digars = push_digar_alt_seq(new_digars, &new_n_digars, &new_m_digar, old_right_digars[i]);
    }
    if (double_check_digar(new_digars, new_n_digars)) {
        // fprintf(stderr, "%s: Error: after updating digar from MSA, invalid digar!\n", tname);
        // print_digar1(new_digars, new_n_digars, stderr);
        // fprintf(stderr, "old digar:\n");
        // print_digar1(digar->digars, digar->n_digar, stderr);
        free_digar1(new_digars, new_n_digars);
    } else {
        free_digar1(digar->digars, digar->n_digar);
        digar->n_digar = new_n_digars; digar->m_digar = new_m_digar; digar->digars = new_digars;
    }
    if (old_left_digars) free(old_left_digars); if (old_right_digars) free(old_right_digars); if (msa_digars) free(msa_digars);
}

void update_digars_from_aln_str(bam_chunk_t *chunk, hts_pos_t noisy_reg_beg, hts_pos_t noisy_reg_end, int *read_id_to_full_covers, int *read_reg_beg, int *read_reg_end,
                                aln_str_t **aln_str, int n_cons, int *clu_n_seqs, int **clu_read_ids) {
    for (int i = 0; i < n_cons; ++i) {
        int n_seqs = clu_n_seqs[i];
        for (int read_i = 0; read_i < n_seqs; ++read_i) {
            int read_id = clu_read_ids[i][read_i]; int read_beg = read_reg_beg[read_id], read_end = read_reg_end[read_id];
            aln_str_t *ref_read_aln_str = LONGCALLD_REF_READ_ALN_STR(aln_str[i], read_i);
            update_digars_from_msa1(chunk->tname, chunk->digars+read_id, ref_read_aln_str->aln_len, ref_read_aln_str->target_aln, ref_read_aln_str->query_aln, read_id_to_full_covers[read_id], noisy_reg_beg, noisy_reg_end, read_beg, read_end);
            if (LONGCALLD_VERBOSE >= 2) fprintf(stderr, "%s digar->qlen: %d\n", bam_get_qname(chunk->reads[read_id]), digar2qlen(chunk->digars+read_id));
        }
    }
}

// return n_cons
// 1. consensu calling; 2. WFA-based MSA

extern int collect_noisy_read_info1(const call_var_opt_t *opt, bam_chunk_t *chunk, int read_id, hts_pos_t reg_beg, hts_pos_t reg_end, int *read_len,
                                    uint8_t **read_seqs, int *fully_covers, int *read_reg_beg, int *read_reg_end);

// extend [win_beg, win_end] so that neither boundary cuts through an indel run of this walk-vs-ref alignment
static int deknot_extend_ref_window(aln_str_t *as, hts_pos_t aln_ref_beg, hts_pos_t *win_beg, hts_pos_t *win_end) {
    int extended = 0;
    hts_pos_t ref_pos = aln_ref_beg;
    int col_beg = -1, col_end = -1;
    for (int col = 0; col < as->aln_len; ++col) {
        if (as->target_aln[col] != 5) {
            if (ref_pos == *win_beg) col_beg = col;
            if (ref_pos == *win_end) { col_end = col; break; }
            ref_pos++;
        }
    }
    if (col_beg == -1 || col_end == -1) return 0;
    hts_pos_t new_beg = *win_beg, new_end = *win_end;
    while (col_beg > 0 && as->query_aln[col_beg] == 5) { // boundary inside a deletion run
        col_beg--;
        if (as->target_aln[col_beg] != 5) new_beg--;
    }
    while (col_end < as->aln_len-1 && as->query_aln[col_end] == 5) {
        col_end++;
        if (as->target_aln[col_end] != 5) new_end++;
    }
    if (new_beg < *win_beg) { *win_beg = new_beg; extended = 1; }
    if (new_end > *win_end) { *win_end = new_end; extended = 1; }
    return extended;
}

// cut the alignment down to the columns covering ref window [win_beg, win_end]
static void deknot_trim_aln(aln_str_t *as, hts_pos_t aln_ref_beg, hts_pos_t win_beg, hts_pos_t win_end) {
    hts_pos_t ref_pos = aln_ref_beg;
    int col_beg = -1, col_end = -1;
    for (int col = 0; col < as->aln_len; ++col) {
        if (as->target_aln[col] != 5) {
            if (ref_pos == win_beg && col_beg == -1) col_beg = col;
            if (ref_pos == win_end) { col_end = col; break; }
            ref_pos++;
        }
    }
    if (col_beg == -1) col_beg = 0;
    if (col_end == -1) col_end = as->aln_len-1;
    int new_len = col_end - col_beg + 1;
    // single block, query_aln points inside it (same convention as wfa_collect_pretty_alignment,
    // downstream cleanup frees target_aln only)
    uint8_t *mem = (uint8_t*)malloc(2 * new_len * sizeof(uint8_t));
    memcpy(mem, as->target_aln + col_beg, new_len);
    memcpy(mem + new_len, as->query_aln + col_beg, new_len);
    free(as->target_aln);
    as->target_aln = mem; as->query_aln = mem + new_len; as->aln_len = new_len;
    as->target_beg = 0; as->target_end = new_len-1;
    as->query_beg = 0; as->query_end = new_len-1;
}

#define LONGCALLD_DEKNOT_SPAN_FLANK 500     // a usable walk aligns through the noisy region with this many bases on both sides
#define LONGCALLD_DEKNOT_MAX_ALN_WALKS 8    // with more than two walks, only the longest this many are aligned
#define LONGCALLD_DEKNOT_IDENT_MAX_ED 5     // two walks are near-identical over the noisy region when their edit distance
#define LONGCALLD_DEKNOT_IDENT_MAX_FRAC 0.005 // is at most max(IDENT_MAX_ED, IDENT_MAX_FRAC * length) with no SV-sized gap

// edlib end-to-end alignment of two sequences: returns the edit distance (-1 on failure) and the longest gap run (bp)
static int deknot_pair_diff(uint8_t *a, int alen, uint8_t *b, int blen, int *max_gap_run) {
    *max_gap_run = 0;
    if (alen <= 0 || blen <= 0) { *max_gap_run = (alen > blen) ? alen : blen; return *max_gap_run; }
    EdlibAlignResult r = edlibAlign((const char*)a, alen, (const char*)b, blen,
                                    edlibNewAlignConfig(-1, EDLIB_MODE_NW, EDLIB_TASK_PATH, NULL, 0));
    if (r.status != EDLIB_STATUS_OK || r.alignment == NULL) { edlibFreeAlignResult(r); return -1; }
    int run = 0;
    for (int i = 0; i < r.alignmentLength; ++i) {
        if (r.alignment[i] == EDLIB_EDOP_INSERT || r.alignment[i] == EDLIB_EDOP_DELETE) {
            run = (i > 0 && r.alignment[i-1] == r.alignment[i]) ? run + 1 : 1;
            if (run > *max_gap_run) *max_gap_run = run;
        }
    }
    int ed = r.editDistance;
    edlibFreeAlignResult(r);
    return ed;
}

// query bases of a walk-vs-ref alignment that fall within ref offsets [beg, end] (insertions inside the range included)
static uint8_t *deknot_aln_sub_query(aln_str_t *as, int beg, int end, int *len) {
    uint8_t *s = (uint8_t*)malloc((as->aln_len > 0 ? as->aln_len : 1) * sizeof(uint8_t));
    *len = 0;
    int ref_pos = 0; // offset of the next ref base
    for (int c = 0; c < as->aln_len; ++c) {
        int t = as->target_aln[c], q = as->query_aln[c];
        if (t != 5) {
            if (ref_pos > end) break;
            if (ref_pos >= beg && q != 5) s[(*len)++] = q;
            ref_pos++;
        } else if (ref_pos > beg && ref_pos <= end && q != 5) s[(*len)++] = q;
    }
    return s;
}

#define LONGCALLD_DEKNOT_PAIR_EDGE_GAP 1000

static int deknot_aln_edge_gap(aln_str_t *as, int min_len, int edge) {
    int c = 0;
    while (c < as->aln_len) {
        int t = as->target_aln[c], q = as->query_aln[c];
        int kind = (t == 5) ? 1 : (q == 5) ? 2 : 0;
        if (kind == 0) { c++; continue; }
        int beg = c;
        while (c < as->aln_len && ((kind == 1 && as->target_aln[c] == 5) || (kind == 2 && as->query_aln[c] == 5 && as->target_aln[c] != 5))) c++;
        int len = c - beg;
        if (len >= min_len && (beg < edge || as->aln_len - c < edge)) return 1;
    }
    return 0;
}

static hts_pos_t deknot_aln_net_change(aln_str_t *as, int beg, int end) {
    hts_pos_t del = 0, ins = 0; int ref_pos = 0;
    for (int c = 0; c < as->aln_len; ++c) {
        int t = as->target_aln[c], q = as->query_aln[c];
        if (t != 5) {
            if (ref_pos > end) break;
            if (ref_pos >= beg && q == 5) del++;
            ref_pos++;
        } else if (ref_pos > beg && ref_pos <= end) ins++;
    }
    return del - ins;
}

// orient and align one DeKnot walk end-to-end to the expanded reference. returns 0 (and an empty *as) when the walk
// aligns poorly in both orientations. on success fills *as, the ref offsets [ref_beg, ref_end] of its aligned bases
// and a score (matches minus mismatches minus gap runs, so an SV allele costs one)
static int deknot_align_walk(const call_var_opt_t *opt, uint8_t *ref_seq, int ref_seq_len, const char *walk, aln_str_t *as,
                             int *ref_beg, int *ref_end, int *score, double cap_div, int *capped) {
    as->target_aln = NULL; as->query_aln = NULL; as->aln_len = 0;
    *ref_beg = -1; *ref_end = -1; *score = 0;
    int wlen = strlen(walk);
    if (wlen == 0) return 0;
    uint8_t *wseq = (uint8_t*)malloc(wlen);
    uint8_t *wseq_rc = (uint8_t*)malloc(wlen);
    for (int j = 0; j < wlen; ++j) wseq[j] = nst_nt4_table[(int)walk[j]];
    for (int j = 0; j < wlen; ++j) wseq_rc[j] = (wseq[wlen-1-j] < 4) ? 3 - wseq[wlen-1-j] : 4;
    // DeKnot walks have arbitrary strand: pick the orientation with more matched bases,
    // and reject walks that align poorly either way
    int fwd_eq = 0, fwd_xid = 0, rev_eq = 0, rev_xid = 0;
    edlib_end2end_aln(ref_seq, ref_seq_len, wseq, wlen, &fwd_eq, &fwd_xid);
    edlib_end2end_aln(ref_seq, ref_seq_len, wseq_rc, wlen, &rev_eq, &rev_xid);
    uint8_t *use_seq = (rev_eq > fwd_eq) ? wseq_rc : wseq;
    int best_eq = (rev_eq > fwd_eq) ? rev_eq : fwd_eq;
    int min_len = (wlen < ref_seq_len) ? wlen : ref_seq_len;
    if (best_eq < (int)(0.55 * min_len)) { free(wseq); free(wseq_rc); return 0; }
    if (capped != NULL) *capped = 0;
    if (cap_div > 0) {
        long long d = wlen > ref_seq_len ? wlen - ref_seq_len : ref_seq_len - wlen;
        long long g1 = opt->gap_open1 + (long long)opt->gap_ext1 * d, g2 = opt->gap_open2 + (long long)opt->gap_ext2 * d;
        long long budget = (d > 0 ? (g1 < g2 ? g1 : g2) : 0) + (long long)(cap_div * min_len * opt->mismatch) + 1000;
        wfa_max_steps = budget < INT_MAX ? (int)budget : INT_MAX;
    }
    wfa_end2end_aln_mem(1, 0, ref_seq, ref_seq_len, use_seq, wlen, opt->gap_aln, opt->mismatch, opt->gap_open1, opt->gap_ext1, opt->gap_open2, opt->gap_ext2,
                        LONGCALLD_WFA_NO_HEURISTIC, LONGCALLD_WFA_AFFINE_2P,
                        NULL, NULL, &as->target_aln, &as->query_aln, &as->aln_len);
    wfa_max_steps = INT_MAX;
    if (capped != NULL) *capped = wfa_last_capped;
    free(wseq); free(wseq_rc);
    if (as->aln_len <= 0) {
        if (as->target_aln != NULL) free(as->target_aln);
        as->target_aln = NULL; as->query_aln = NULL; as->aln_len = 0;
        return 0;
    }
    as->target_beg = 0; as->target_end = as->aln_len-1;
    as->query_beg = 0; as->query_end = as->aln_len-1;
    int ref_pos = 0, eq = 0, xid = 0, gaps = 0, prev_gap = 0;
    for (int c = 0; c < as->aln_len; ++c) {
        int t = as->target_aln[c], q = as->query_aln[c];
        int cur_gap = (t == 5) ? 1 : (q == 5) ? 2 : 0;
        if (cur_gap == 0) {
            if (*ref_beg < 0) *ref_beg = ref_pos;
            *ref_end = ref_pos;
            if (t == q) eq++; else xid++;
        } else if (cur_gap != prev_gap) gaps++;
        prev_gap = cur_gap;
        if (t != 5) ref_pos++;
    }
    *score = eq - xid - gaps;
    return 1;
}

#define LONGCALLD_DEKNOT_PAIR_MIN_SV 5000 // a breakpoint pair predicts an indel of at least this, or half the breakpoint distance
#define LONGCALLD_DEKNOT_PAIR_WFA_CAP_DIV 0.02

// longest indel run of at least min_sv bp in a walk-vs-ref alignment (ref coordinate of column 0: aln_ref_beg) whose
// reference span overlaps [gap_beg, gap_end]; returns its length (0: none) and span (a deletion: the deleted
// reference bases; an insertion: the two reference bases around it)
#define LONGCALLD_DEKNOT_PAIR_SV_ANCHOR 200
static int deknot_walk_sv_span2(aln_str_t *as, hts_pos_t aln_ref_beg, hts_pos_t gap_beg, hts_pos_t gap_end, hts_pos_t min_sv,
                                hts_pos_t *sv_beg, hts_pos_t *sv_end, int *q_at) {
    hts_pos_t ref_pos = aln_ref_beg; // reference position of the next reference base
    int best = 0, c = 0, qpos = 0;
    int n_aligned = 0, aligned_before = 0;
    for (int k = 0; k < as->aln_len; ++k) if (as->target_aln[k] != 5 && as->query_aln[k] != 5) n_aligned++;
    while (c < as->aln_len) {
        int t = as->target_aln[c], q = as->query_aln[c];
        if (t != 5 && q == 5) { // deletion run
            hts_pos_t rb = ref_pos; int len = 0, q0 = qpos;
            while (c < as->aln_len && as->target_aln[c] != 5 && as->query_aln[c] == 5) { len++; ref_pos++; c++; }
            hts_pos_t re = ref_pos - 1;
            int anchored = (aligned_before >= LONGCALLD_DEKNOT_PAIR_SV_ANCHOR && n_aligned - aligned_before >= LONGCALLD_DEKNOT_PAIR_SV_ANCHOR);
            if (anchored && len >= min_sv && len > best && rb <= gap_end && re >= gap_beg) { best = len; *sv_beg = rb; *sv_end = re; if (q_at) *q_at = q0; }
        } else if (t == 5) { // insertion run
            int len = 0, q0 = qpos;
            while (c < as->aln_len && as->target_aln[c] == 5) { len++; qpos++; c++; }
            hts_pos_t rb = ref_pos - 1, re = ref_pos;
            int anchored = (aligned_before >= LONGCALLD_DEKNOT_PAIR_SV_ANCHOR && n_aligned - aligned_before >= LONGCALLD_DEKNOT_PAIR_SV_ANCHOR);
            if (anchored && len >= min_sv && len > best && rb <= gap_end && re >= gap_beg) { best = len; *sv_beg = rb; *sv_end = re; if (q_at) *q_at = q0; }
        } else { ref_pos++; if (q != 5) { qpos++; aligned_before++; } c++; }
    }
    return best;
}
static int deknot_walk_sv_span(aln_str_t *as, hts_pos_t aln_ref_beg, hts_pos_t gap_beg, hts_pos_t gap_end, hts_pos_t min_sv,
                               hts_pos_t *sv_beg, hts_pos_t *sv_end) {
    return deknot_walk_sv_span2(as, aln_ref_beg, gap_beg, gap_end, min_sv, sv_beg, sv_end, NULL);
}

int32_t deknot_all_depth_median(bam_chunk_t *chunk) {
    if (chunk->all_depth == NULL || chunk->n_all_depth <= 0) return 0;
    int n = chunk->n_all_depth;
    int32_t *v = (int32_t*)malloc(n * sizeof(int32_t));
    memcpy(v, chunk->all_depth, n * sizeof(int32_t));
    for (int i = 1; i < n; ++i) { int32_t x = v[i]; int j = i; while (j > 0 && v[j-1] > x) { v[j] = v[j-1]; j--; } v[j] = x; }
    int32_t m = v[n / 2];
    free(v);
    return m;
}

#define LONGCALLD_DEKNOT_UNIQUE_SEG 2000
#define LONGCALLD_DEKNOT_UNIQUE_SEG_DIV 0.10
#define LONGCALLD_DEKNOT_UNIQUE_SEARCH 200000

static int deknot_ref_seg_unique(bam_chunk_t *chunk, hts_pos_t seg_beg) {
    int L = LONGCALLD_DEKNOT_UNIQUE_SEG;
    if (seg_beg < chunk->ref_beg || seg_beg + L - 1 > chunk->ref_end) return 0;
    hts_pos_t lo = seg_beg - LONGCALLD_DEKNOT_UNIQUE_SEARCH, hi = seg_beg + L - 1 + LONGCALLD_DEKNOT_UNIQUE_SEARCH;
    if (lo < chunk->ref_beg) lo = chunk->ref_beg;
    if (hi > chunk->ref_end) hi = chunk->ref_end;
    int tlen = (int)(hi - lo + 1);
    char *seg = (char*)malloc(L), *tgt = (char*)malloc(tlen);
    int n_n = 0;
    for (int i = 0; i < L; ++i) { seg[i] = toupper(chunk->ref_seq[seg_beg + i - chunk->ref_beg]); if (seg[i] == 'N') n_n++; }
    for (int i = 0; i < tlen; ++i) tgt[i] = toupper(chunk->ref_seq[lo + i - chunk->ref_beg]);
    for (hts_pos_t p = seg_beg; p < seg_beg + L; ++p) tgt[p - lo] = 'N';
    int uniq = 0;
    if (n_n * 10 < L) {
        EdlibAlignResult r = edlibAlign(seg, L, tgt, tlen, edlibNewAlignConfig((int)(LONGCALLD_DEKNOT_UNIQUE_SEG_DIV * L), EDLIB_MODE_HW, EDLIB_TASK_DISTANCE, NULL, 0));
        uniq = (r.status == EDLIB_STATUS_OK && r.editDistance < 0);
        edlibFreeAlignResult(r);
    }
    free(seg); free(tgt);
    return uniq;
}

hts_pos_t deknot_find_unique_flank(bam_chunk_t *chunk, hts_pos_t from, int dir, hts_pos_t limit, int32_t median) {
    if (chunk->hq_depth == NULL || median <= 0) return -1;
    double thr = LONGCALLD_DEKNOT_UNIQUE_FRAC * median;
    int need = LONGCALLD_DEKNOT_UNIQUE_LEN / LONGCALLD_DEPTH_BIN;
    int b = (int)((from - chunk->reg_beg) / LONGCALLD_DEPTH_BIN);
    int lb = (int)((limit - chunk->reg_beg) / LONGCALLD_DEPTH_BIN);
    int run = 0, run_inner = -1;
    for (; dir < 0 ? b >= lb : b <= lb; b += dir) {
        if (b < 0 || b >= chunk->n_all_depth) break;
        if (chunk->hq_depth[b] >= thr) {
            if (run == 0) run_inner = b;
            if (++run >= need) {
                hts_pos_t edge = chunk->reg_beg + (hts_pos_t)run_inner * LONGCALLD_DEPTH_BIN;
                hts_pos_t inner = dir < 0 ? edge + LONGCALLD_DEPTH_BIN - 1 : edge;
                hts_pos_t seg_beg = dir < 0 ? inner - LONGCALLD_DEKNOT_UNIQUE_SEG + 1 : inner;
                if (deknot_ref_seg_unique(chunk, seg_beg)) return inner;
                run = 0;
            }
        } else run = 0;
    }
    return -1;
}

static int deknot_junction_reads(const call_var_opt_t *opt, aln_str_t *as, int q_at, char **read_strs, int n_reads) {
    int qlen = 0;
    for (int c = 0; c < as->aln_len; ++c) if (as->query_aln[c] != 5) qlen++;
    int f = LONGCALLD_DEKNOT_JUNC_FLANK;
    if (q_at - f < 0 || q_at + f > qlen) return 0;
    char *seg = (char*)malloc(2 * f);
    int qi = 0, k = 0;
    for (int c = 0; c < as->aln_len && k < 2 * f; ++c) {
        if (as->query_aln[c] == 5) continue;
        if (qi >= q_at - f && qi < q_at + f) seg[k++] = "ACGTN"[as->query_aln[c] < 5 ? as->query_aln[c] : 4];
        qi++;
    }
    int max_ed = (int)((opt->is_ont ? 0.10 : 0.03) * 2 * f), n = 0;
    for (int i = 0; i < n_reads; ++i) {
        int rlen = strlen(read_strs[i]);
        if (rlen < 2 * f) continue;
        EdlibAlignResult r = edlibAlign(seg, 2 * f, read_strs[i], rlen, edlibNewAlignConfig(max_ed, EDLIB_MODE_HW, EDLIB_TASK_DISTANCE, NULL, 0));
        if (r.status == EDLIB_STATUS_OK && r.editDistance >= 0) n++;
        edlibFreeAlignResult(r);
    }
    free(seg);
    return n;
}

static int deknot_merge_large_dels(aln_str_t *as, int min_len) {
    int c0 = -1, c1 = -1, n_runs = 0, c = 0;
    while (c < as->aln_len) {
        if (as->target_aln[c] != 5 && as->query_aln[c] == 5) {
            int b = c;
            while (c < as->aln_len && as->target_aln[c] != 5 && as->query_aln[c] == 5) c++;
            if (c - b >= min_len) { if (c0 < 0) c0 = b; c1 = c - 1; n_runs++; }
        } else c++;
    }
    if (n_runs < 2) return 0;
    int lr = 0, lq = 0;
    uint8_t *R = (uint8_t*)malloc(c1 - c0 + 1), *Q = (uint8_t*)malloc(c1 - c0 + 1);
    for (c = c0; c <= c1; ++c) {
        if (as->target_aln[c] != 5) R[lr++] = as->target_aln[c];
        if (as->query_aln[c] != 5) Q[lq++] = as->query_aln[c];
    }
    if (lq >= lr) { free(R); free(Q); return 0; }
    int del = lr - lq, left = 1;
    EdlibAlignResult best;
    if (lq > 0) {
        EdlibAlignResult rl = edlibAlign((char*)Q, lq, (char*)R, lq, edlibNewAlignConfig(-1, EDLIB_MODE_NW, EDLIB_TASK_PATH, NULL, 0));
        EdlibAlignResult rr = edlibAlign((char*)Q, lq, (char*)R + del, lq, edlibNewAlignConfig(-1, EDLIB_MODE_NW, EDLIB_TASK_PATH, NULL, 0));
        if (rl.status != EDLIB_STATUS_OK || rr.status != EDLIB_STATUS_OK) { edlibFreeAlignResult(rl); edlibFreeAlignResult(rr); free(R); free(Q); return 0; }
        if (rr.editDistance < rl.editDistance) { left = 0; best = rr; edlibFreeAlignResult(rl); }
        else { best = rl; edlibFreeAlignResult(rr); }
    } else { best.alignment = NULL; best.alignmentLength = 0; }
    int mid = best.alignmentLength + del;
    int new_len = c0 + mid + (as->aln_len - c1 - 1);
    uint8_t *mem = (uint8_t*)malloc(2 * new_len);
    uint8_t *nt = mem, *nq = mem + new_len;
    memcpy(nt, as->target_aln, c0); memcpy(nq, as->query_aln, c0);
    int o = c0, ti = left ? 0 : del, qi = 0;
    if (!left) for (int i = 0; i < del; ++i, ++o) { nt[o] = R[i]; nq[o] = 5; }
    for (int i = 0; i < best.alignmentLength; ++i, ++o) {
        unsigned char op = best.alignment[i];
        if (op == EDLIB_EDOP_INSERT) { nt[o] = 5; nq[o] = Q[qi++]; }
        else if (op == EDLIB_EDOP_DELETE) { nt[o] = R[ti++]; nq[o] = 5; }
        else { nt[o] = R[ti++]; nq[o] = Q[qi++]; }
    }
    if (left) for (int i = lq; i < lr; ++i, ++o) { nt[o] = R[i]; nq[o] = 5; }
    memcpy(nt + o, as->target_aln + c1 + 1, as->aln_len - c1 - 1);
    memcpy(nq + o, as->query_aln + c1 + 1, as->aln_len - c1 - 1);
    if (lq > 0) edlibFreeAlignResult(best);
    free(R); free(Q);
    free(as->target_aln);
    as->target_aln = nt; as->query_aln = nq; as->aln_len = new_len;
    as->target_beg = 0; as->target_end = new_len - 1; as->query_beg = 0; as->query_end = new_len - 1;
    return 1;
}

static int deknot_large_del_union(aln_str_t *as, hts_pos_t aln_ref_beg, int min_len, hts_pos_t *b, hts_pos_t *e) {
    hts_pos_t ref_pos = aln_ref_beg; int n = 0, c = 0;
    while (c < as->aln_len) {
        if (as->target_aln[c] != 5 && as->query_aln[c] == 5) {
            hts_pos_t rb = ref_pos; int len = 0;
            while (c < as->aln_len && as->target_aln[c] != 5 && as->query_aln[c] == 5) { len++; ref_pos++; c++; }
            if (len >= min_len) { if (n == 0 || rb < *b) *b = rb; if (n == 0 || ref_pos - 1 > *e) *e = ref_pos - 1; n++; }
        } else { if (as->target_aln[c] != 5) ref_pos++; c++; }
    }
    return n;
}

// local assembly rescue over an expanded window: full-length read subsequences are assembled with DeKnot,
// the phased walks are aligned to the expanded reference, and the alignments are trimmed back to the
// noisy region (extended over any indel run that crosses its boundary, so SVs survive intact).
// a result with more than two walks (several connected components) is reduced to the walks that span the
// noisy region, and from those to a pair differing by a large indel or the single best-aligning walk;
// a two-walk result whose walks are near-identical over the noisy region is collapsed to one consensus.
// commits results into clu_n_seqs/clu_read_ids/aln_strs only when it improves on n_cons_existing.
// pair_beg/pair_end >= 0: the partner breakpoint of the noisy region (both breakpoints of one large SV are known);
// the window, the span check and the trimmed region then cover both breakpoints.
static int deknot_lowmq_window_seq(bam1_t *b, hts_pos_t beg, hts_pos_t end, char **out) {
    uint32_t *cigar = bam_get_cigar(b); int n_cigar = b->core.n_cigar, qlen = b->core.l_qseq;
    hts_pos_t ref = b->core.pos + 1; int q = 0, qs = -1, qe = -1;
    for (int i = 0; i < n_cigar; ++i) {
        int op = bam_cigar_op(cigar[i]), len = bam_cigar_oplen(cigar[i]);
        if (op == BAM_CMATCH || op == BAM_CEQUAL || op == BAM_CDIFF) {
            if (qs < 0 && beg >= ref && beg < ref + len) qs = q + (int)(beg - ref);
            if (qe < 0 && end >= ref && end < ref + len) qe = q + (int)(end - ref);
            ref += len; q += len;
        } else if (op == BAM_CINS || op == BAM_CSOFT_CLIP) {
            q += len;
        } else if (op == BAM_CDEL || op == BAM_CREF_SKIP) {
            if (qs < 0 && beg >= ref && beg < ref + len) qs = q;
            if (qe < 0 && end >= ref && end < ref + len) qe = q - 1;
            ref += len;
        }
    }
    if (b->core.pos + 1 > end || ref - 1 < beg) return 0;
    if (qs < 0) qs = 0;
    if (qe < 0) qe = qlen - 1;
    if (qe < qs || qlen <= 0) return 0;
    uint8_t *s = bam_get_seq(b);
    int n = qe - qs + 1;
    *out = (char*)malloc(n + 1);
    for (int i = 0; i < n; ++i) (*out)[i] = seq_nt16_str[bam_seqi(s, qs + i)];
    (*out)[n] = '\0';
    return n;
}

static int deknot_sv_rescue_core(const call_var_opt_t *opt, bam_chunk_t *chunk, hts_pos_t noisy_reg_beg, hts_pos_t noisy_reg_end,
                     hts_pos_t pair_beg, hts_pos_t pair_end,
                     int n_noisy_reads, int *noisy_reads, hts_pos_t flank_len, int abort_on_edge, int *hit_window_edge, int pass,
                     int n_cons_existing, int *clu_n_seqs, int **clu_read_ids, aln_str_t **aln_strs, hts_pos_t *var_reg_beg,
                     hts_pos_t *sv_beg, hts_pos_t *sv_end, hts_pos_t win_lo, hts_pos_t win_hi, int *junc_reads) {
    // region whose variants are wanted: the noisy region, or the span over both breakpoints
    int paired = (pair_beg >= 0 && pair_end >= pair_beg);
    hts_pos_t reg_beg = noisy_reg_beg, reg_end = noisy_reg_end;
    if (paired) {
        if (pair_beg < reg_beg) reg_beg = pair_beg;
        if (pair_end > reg_end) reg_end = pair_end;
    }
    hts_pos_t max_win = pass == DK_PASS_DEPTH ? LONGCALLD_DEKNOT_DEPTH_MAX_WIN : paired ? LONGCALLD_DEKNOT_PAIR_MAX_WIN : LONGCALLD_DEKNOT_MAX_WIN;
    if (junc_reads != NULL) *junc_reads = -1;
    int max_asm_reads = paired ? LONGCALLD_DEKNOT_PAIR_MAX_ASM_READS : LONGCALLD_DEKNOT_MAX_ASM_READS;
    // a region that already has a consensus is only worth re-assembling when it shows a
    // breakpoint signature: multiple reads that fail to fully span it (clips/gaps at the edges)
    if (n_cons_existing >= 1) {
        int n_nonfull = 0;
        for (int k = 0; k < n_noisy_reads; ++k) {
            int rlen = 0, rcover = 0, rb = 0, re = 0; uint8_t *rseq = NULL;
            collect_noisy_read_info1(opt, chunk, noisy_reads[k], noisy_reg_beg, noisy_reg_end, &rlen, &rseq, &rcover, &rb, &re);
            free(rseq);
            if (!LONGCALLD_NOISY_IS_BOTH_COVER(rcover)) n_nonfull++;
        }
        if (n_nonfull < 2) { dk_stats_count(&dk_stats.n_skip_no_breakpoint); return -1; }
    }
    *hit_window_edge = 0;
    hts_pos_t exp_beg = (win_lo >= 0 && win_lo < reg_beg ? win_lo : reg_beg) - flank_len;
    hts_pos_t exp_end = (win_hi > reg_end ? win_hi : reg_end) + flank_len;
    if (exp_beg < 1) exp_beg = 1;
    uint8_t *ref_seq = NULL;
    int ref_seq_len = collect_reg_ref_bseq(chunk, &exp_beg, &exp_end, &ref_seq);
    if (ref_seq_len <= 0 || ref_seq_len > max_win) {
        dk_stats_count(&dk_stats.n_skip_window);
        if (paired) dk_stats_count(&dk_stats.n_pair_pre_reject);
        free(ref_seq);
        return paired ? -1 : n_cons_existing;
    }

    int n_reads = 0;
    int *read_ids = (int*)malloc(chunk->n_reads * sizeof(int));
    for (int i = 0; i < chunk->n_reads; ++i) {
        if (chunk->is_skipped[i]) continue;
        digar_t *d = chunk->digars + i;
        if (d->beg <= exp_end && d->end >= exp_beg) read_ids[n_reads++] = i;
    }
    if (n_reads < opt->min_dp || n_reads > opt->max_noisy_reg_cov) {
        dk_stats_count(&dk_stats.n_skip_depth);
        if (paired) dk_stats_count(&dk_stats.n_pair_pre_reject);
        free(read_ids); free(ref_seq); return -1;
    }
    // assembly does not need full depth: keep every noisy-region read, subsample the rest evenly
    int subsampled = 0;
    if (n_reads > max_asm_reads) {
        subsampled = 1;
        int *kept = (int*)malloc(n_reads * sizeof(int));
        int n_kept = 0;
        for (int i = 0; i < n_reads; ++i) {
            for (int k = 0; k < n_noisy_reads; ++k) {
                if (noisy_reads[k] == read_ids[i]) { kept[n_kept++] = read_ids[i]; read_ids[i] = -1; break; }
            }
        }
        int n_rest = 0;
        for (int i = 0; i < n_reads; ++i) if (read_ids[i] >= 0) n_rest++;
        int want = max_asm_reads - n_kept;
        if (want > 0 && n_rest > 0) {
            double step = (want >= n_rest) ? 1.0 : (double)n_rest / (double)want;
            double next = 0.0; int seen = 0;
            for (int i = 0; i < n_reads && n_kept < max_asm_reads; ++i) {
                if (read_ids[i] < 0) continue;
                if (seen >= (int)next) { kept[n_kept++] = read_ids[i]; next += step; }
                seen++;
            }
        }
        free(read_ids);
        read_ids = kept; n_reads = n_kept;
    }

    int n_lowmq = 0;
    for (int i = 0; i < chunk->n_lowmq_reads; ++i) {
        bam1_t *b = chunk->lowmq_reads[i];
        if (b->core.pos + 1 <= exp_end && bam_endpos(b) >= exp_beg) n_lowmq++;
    }
    char **read_strs = (char**)malloc((n_reads + n_lowmq) * sizeof(char*));
    char **read_names = (char**)malloc((n_reads + n_lowmq) * sizeof(char*));
    for (int i = 0; i < n_reads; ++i) {
        int len = 0, cover = 0, rb = 0, re = 0; uint8_t *bseq = NULL;
        collect_noisy_read_info1(opt, chunk, read_ids[i], exp_beg, exp_end, &len, &bseq, &cover, &rb, &re);
        read_strs[i] = (char*)malloc(len + 1);
        for (int j = 0; j < len; ++j) read_strs[i][j] = "ACGTN"[bseq[j]];
        read_strs[i][len] = '\0';
        free(bseq);
        read_names[i] = (char*)malloc(16);
        snprintf(read_names[i], 16, "r%d", i);
    }
    if (n_lowmq > 0) {
        int n0 = n_reads;
        for (int i = 0; i < chunk->n_lowmq_reads; ++i) {
            bam1_t *b = chunk->lowmq_reads[i];
            if (b->core.pos + 1 > exp_end || bam_endpos(b) < exp_beg) continue;
            char *s = NULL;
            if (deknot_lowmq_window_seq(b, exp_beg, exp_end, &s) <= 0) continue;
            read_strs[n_reads] = s;
            read_names[n_reads] = (char*)malloc(16);
            snprintf(read_names[n_reads], 16, "q%d", n_reads);
            n_reads++;
        }
        if (LONGCALLD_VERBOSE >= 1)
            fprintf(stderr, "DeKnotLowMQ %s:%" PRIi64 "-%" PRIi64 " window %" PRIi64 "-%" PRIi64 " added %d low-MAPQ reads to %d\n",
                    chunk->tname, noisy_reg_beg, noisy_reg_end, exp_beg, exp_end, n_reads - n0, n0);
    }
    dk_opt_t dk_opt; dk_opt_init(&dk_opt);
    dk_opt.check_strand = (opt->is_ont != 0);
    double asm_t0 = realtime();
    dk_walks_t *walks = dk_assemble_walks(&dk_opt, read_strs, read_names, n_reads);
    double asm_secs = realtime() - asm_t0;
    int n_walks = (walks == NULL) ? 0 : walks->n_walks;
    int walks_null = (walks == NULL), demoted = 0;
    int sel = DK_SEL_NONE, collapsed = 0;
    int outcome;

    int rescue_n = 0;
    if (n_walks >= 1 && n_walks <= 2) rescue_n = n_walks;
    // with one or two walks the count is final; with more, the selection below decides it
    if (n_walks == 0 || (n_walks <= 2 && rescue_n <= n_cons_existing)) {
        outcome = (n_walks == 0) ? DK_OUT_NO_WALKS : DK_OUT_NO_GAIN;
        dk_stats_add_attempt(chunk, noisy_reg_beg, noisy_reg_end, exp_beg, exp_end, n_reads, n_noisy_reads, subsampled, n_walks, walks_null, demoted, sel, collapsed, asm_secs, pass, outcome);
        if (walks != NULL) dk_walks_free(walks);
        for (int i = 0; i < n_reads; ++i) { free(read_strs[i]); free(read_names[i]); }
        free(read_strs); free(read_names); free(read_ids); free(ref_seq);
        return -1;
    }

    // a usable walk aligns through the noisy region with LONGCALLD_DEKNOT_SPAN_FLANK bases on both sides
    // (ref offsets within the expanded window)
    hts_pos_t span_beg = reg_beg - LONGCALLD_DEKNOT_SPAN_FLANK, span_end = reg_end + LONGCALLD_DEKNOT_SPAN_FLANK;
    if (span_beg < exp_beg) span_beg = exp_beg;
    if (span_end > exp_end) span_end = exp_end;
    int span_beg_off = (int)(span_beg - exp_beg), span_end_off = (int)(span_end - exp_beg);

    // paired pass: the walks have to deliver the SV the breakpoint pair predicts, an indel of at least
    // max(5 kb, half the breakpoint distance) whose reference span overlaps the interval between the two breakpoints
    hts_pos_t gap_beg = 0, gap_end = 0, min_sv = 0;
    if (paired) {
        if (pair_beg > noisy_reg_end) { gap_beg = noisy_reg_end; gap_end = pair_beg; }
        else { gap_beg = pair_end; gap_end = noisy_reg_beg; }
        min_sv = (gap_end - gap_beg) / 2;
        if (min_sv < LONGCALLD_DEKNOT_PAIR_MIN_SV) min_sv = LONGCALLD_DEKNOT_PAIR_MIN_SV;
    }
    // walks to align: all of them for one or two walks; for more, the longest few that are long enough to span
    int n_aln = 0;
    int *aln_order = (int*)malloc(n_walks * sizeof(int));
    int *wlens = (int*)malloc(n_walks * sizeof(int));
    for (int i = 0; i < n_walks; ++i) wlens[i] = strlen(walks->walks[i]);
    if (n_walks <= 2) {
        for (int i = 0; i < n_walks; ++i) aln_order[n_aln++] = i;
    } else {
        // a walk shorter than the span cannot cross it, except over a breakpoint pair where a deletion walk is
        // shorter than the reference span by the deleted length: there only the two flanks are required
        int min_wlen = paired ? 2 * LONGCALLD_DEKNOT_SPAN_FLANK : span_end_off - span_beg_off + 1;
        for (int i = 0; i < n_walks; ++i) { // insertion sort by length, longest first
            if (wlens[i] < min_wlen) continue;
            int k = n_aln++;
            while (k > 0 && wlens[aln_order[k-1]] < wlens[i]) { aln_order[k] = aln_order[k-1]; k--; }
            aln_order[k] = i;
        }
        if (n_aln > LONGCALLD_DEKNOT_MAX_ALN_WALKS) n_aln = LONGCALLD_DEKNOT_MAX_ALN_WALKS;
    }
    aln_str_t *cand_alns = (aln_str_t*)calloc(n_aln > 0 ? n_aln : 1, sizeof(aln_str_t));
    int *cand_score = (int*)malloc((n_aln > 0 ? n_aln : 1) * sizeof(int));
    int n_cand = 0, aln_ok = 1;
    double aln_t0 = realtime();
    for (int k = 0; k < n_aln; ++k) {
        int rb = -1, re = -1, sc = 0;
        int capped = 0;
        int ok = deknot_align_walk(opt, ref_seq, ref_seq_len, walks->walks[aln_order[k]], cand_alns + n_cand, &rb, &re, &sc,
                                   paired ? LONGCALLD_DEKNOT_PAIR_WFA_CAP_DIV : 0.0, &capped);
        if (capped) dk_stats_count(&dk_stats.n_wfa_capped);
        if (LONGCALLD_VERBOSE >= 2 || (paired && LONGCALLD_VERBOSE >= 1))
            fprintf(stderr, "DeKnotWalk %s:%" PRIi64 "-%" PRIi64 " walk=%d len=%d aligned=%d%s ref=%" PRIi64 "-%" PRIi64 " score=%d span=%d\n",
                    chunk->tname, noisy_reg_beg, noisy_reg_end, aln_order[k], wlens[aln_order[k]], ok, capped ? " (WFA score cap reached)" : "",
                    ok ? exp_beg + rb : (hts_pos_t)-1, ok ? exp_beg + re : (hts_pos_t)-1, sc, ok && rb <= span_beg_off && re >= span_end_off);
        if (!ok) {
            if (n_walks <= 2 && !paired) { aln_ok = 0; break; }
            continue;
        }
        if (ok && n_walks > 2 && paired && deknot_aln_edge_gap(cand_alns + n_cand, (int)min_sv, LONGCALLD_DEKNOT_PAIR_EDGE_GAP)) {
            if (LONGCALLD_VERBOSE >= 2) fprintf(stderr, "DeKnotWalk %s:%" PRIi64 "-%" PRIi64 " walk=%d dropped: SV run touches the window edge\n", chunk->tname, noisy_reg_beg, noisy_reg_end, aln_order[k]);
            free(cand_alns[n_cand].target_aln);
            cand_alns[n_cand].target_aln = NULL; cand_alns[n_cand].query_aln = NULL; cand_alns[n_cand].aln_len = 0;
            continue;
        }
        if (n_walks > 2 && (rb > span_beg_off || re < span_end_off)) { // no aligned bases through both flanks
            free(cand_alns[n_cand].target_aln);
            cand_alns[n_cand].target_aln = NULL; cand_alns[n_cand].query_aln = NULL; cand_alns[n_cand].aln_len = 0;
            continue;
        }
        cand_score[n_cand++] = sc;
    }
    free(aln_order); free(wlens);
    dk_aln_secs = realtime() - aln_t0;

    aln_str_t walk_alns[2]; int walk_score[2] = {0, 0};
    for (int i = 0; i < 2; ++i) { walk_alns[i].target_aln = NULL; walk_alns[i].query_aln = NULL; walk_alns[i].aln_len = 0; }
    if (aln_ok && n_walks > 2) {
        // choose at most two spanning walks: the best-scoring pair that differs by an SV-sized indel over the
        // noisy region, otherwise the single best-aligning walk
        int pick0 = -1, pick1 = -1;
        if (n_cand >= 2) {
            uint8_t **sub = (uint8_t**)malloc(n_cand * sizeof(uint8_t*));
            int *sub_len = (int*)malloc(n_cand * sizeof(int));
            hts_pos_t *net = (hts_pos_t*)malloc(n_cand * sizeof(hts_pos_t));
            for (int i = 0; i < n_cand; ++i) {
                sub[i] = deknot_aln_sub_query(cand_alns + i, span_beg_off, span_end_off, sub_len + i);
                net[i] = deknot_aln_net_change(cand_alns + i, span_beg_off, span_end_off);
            }
            int best = INT32_MIN;
            for (int i = 0; i < n_cand; ++i) {
                for (int j = i+1; j < n_cand; ++j) {
                    int gap = 0, ed = deknot_pair_diff(sub[i], sub_len[i], sub[j], sub_len[j], &gap);
                    hts_pos_t dnet = net[i] > net[j] ? net[i] - net[j] : net[j] - net[i];
                    int differ = (ed >= 0 && gap >= opt->min_sv_len) || (paired && dnet >= opt->min_sv_len);
                    if (LONGCALLD_VERBOSE >= 2)
                        fprintf(stderr, "DeKnotPairSel %s:%" PRIi64 "-%" PRIi64 " cand=%d,%d sub_len=%d,%d ed=%d max_gap=%d net=%" PRIi64 ",%" PRIi64 " scores=%d,%d differ=%d\n",
                                chunk->tname, noisy_reg_beg, noisy_reg_end, i, j, sub_len[i], sub_len[j], ed, gap, net[i], net[j], cand_score[i], cand_score[j], differ);
                    if (!differ) continue;
                    if (cand_score[i] + cand_score[j] > best) { best = cand_score[i] + cand_score[j]; pick0 = i; pick1 = j; }
                }
            }
            for (int i = 0; i < n_cand; ++i) free(sub[i]);
            free(sub); free(sub_len); free(net);
        }
        if (pick0 < 0) {
            for (int i = 0; i < n_cand; ++i) if (pick0 < 0 || cand_score[i] > cand_score[pick0]) pick0 = i;
        }
        for (int i = 0; i < n_cand; ++i) {
            if (i == pick0) { walk_alns[0] = cand_alns[i]; walk_score[0] = cand_score[i]; }
            else if (i == pick1) { walk_alns[1] = cand_alns[i]; walk_score[1] = cand_score[i]; }
            else free(cand_alns[i].target_aln);
        }
        rescue_n = (pick1 >= 0) ? 2 : (pick0 >= 0) ? 1 : 0;
        sel = (rescue_n == 2) ? DK_SEL_PAIR : (rescue_n == 1) ? DK_SEL_SINGLE : DK_SEL_NONE;
    } else if (aln_ok) {
        for (int i = 0; i < n_cand; ++i) { walk_alns[i] = cand_alns[i]; walk_score[i] = cand_score[i]; }
        if (paired) rescue_n = n_cand;
    } else {
        for (int i = 0; i < n_cand; ++i) free(cand_alns[i].target_aln);
    }
    free(cand_alns); free(cand_score);
    if (aln_ok && n_walks > 2 && rescue_n <= n_cons_existing) {
        outcome = (rescue_n == 0) ? DK_OUT_MANY_WALKS : DK_OUT_NO_GAIN;
        for (int i = 0; i < rescue_n; ++i) free(walk_alns[i].target_aln);
        dk_stats_add_attempt(chunk, noisy_reg_beg, noisy_reg_end, exp_beg, exp_end, n_reads, n_noisy_reads, subsampled, n_walks, walks_null, demoted, sel, collapsed, asm_secs, pass, outcome);
        dk_walks_free(walks);
        for (int i = 0; i < n_reads; ++i) { free(read_strs[i]); free(read_names[i]); }
        free(read_strs); free(read_names); free(read_ids); free(ref_seq);
        return -1;
    }
    hts_pos_t pred_sv_beg = -1, pred_sv_end = -1; // reference span of the predicted SV over the selected walks
    if (aln_ok && paired) {
        int found = 0, junc_ok = 1, junc_min = -1, keep[2] = {1, 1}, n_dropped = 0, n_sv = 0;
        int32_t median = deknot_all_depth_median(chunk);
        for (int i = 0; i < rescue_n; ++i) {
            hts_pos_t sb = -1, se = -1; int q_at = -1;
            int len = deknot_walk_sv_span2(walk_alns + i, exp_beg, gap_beg, gap_end, min_sv, &sb, &se, &q_at);
            if (len > 0) {
                n_sv++;
                int nj = deknot_junction_reads(opt, walk_alns + i, q_at, read_strs, n_reads);
                hts_pos_t lu = deknot_find_unique_flank(chunk, sb - 1, -1, exp_beg, median);
                hts_pos_t ru = deknot_find_unique_flank(chunk, se + 1, 1, exp_end, median);
                int anchored = (lu >= 0 || ru >= 0 || median <= 0);
                int pass_gate = (nj >= LONGCALLD_DEKNOT_JUNC_MIN_READS && anchored);
                if (LONGCALLD_VERBOSE >= 1)
                    fprintf(stderr, "DeKnotJunction %s:%" PRIi64 "-%" PRIi64 " pass=%s walk=%d sv=%" PRIi64 "-%" PRIi64 " (%d) junction_reads=%d/%d unique_flank_left=%" PRIi64 " right=%" PRIi64 " median=%d%s\n",
                            chunk->tname, noisy_reg_beg, noisy_reg_end, dk_pass_names[pass], i, sb, se, len, nj, n_reads, lu, ru, median,
                            pass_gate ? "" : " DROP");
                if (!pass_gate) { keep[i] = 0; n_dropped++; continue; }
                if (junc_min < 0 || nj < junc_min) junc_min = nj;
            }
            if (LONGCALLD_VERBOSE >= 1)
                fprintf(stderr, "DeKnotPairGate %s:%" PRIi64 "-%" PRIi64 " walk=%d predicted>=%" PRIi64 " in %" PRIi64 "-%" PRIi64 " found=%d span=%" PRIi64 "-%" PRIi64 "\n",
                        chunk->tname, noisy_reg_beg, noisy_reg_end, i, min_sv, gap_beg, gap_end, len, sb, se);
            if (len > 0) {
                found = 1;
                if (pred_sv_beg < 0 || sb < pred_sv_beg) pred_sv_beg = sb;
                if (se > pred_sv_end) pred_sv_end = se;
            }
        }
        if (n_dropped > 0) {
            int n_kept = 0;
            for (int i = 0; i < rescue_n; ++i) {
                if (!keep[i]) { free(walk_alns[i].target_aln); walk_alns[i].target_aln = NULL; continue; }
                if (n_kept != i) { walk_alns[n_kept] = walk_alns[i]; walk_score[n_kept] = walk_score[i]; walk_alns[i].target_aln = NULL; }
                n_kept++;
            }
            rescue_n = n_kept;
        }
        if (!found && n_sv > 0) junc_ok = 0;
        if (junc_reads != NULL) *junc_reads = junc_min;
        if (!found || !junc_ok) {
            for (int i = 0; i < rescue_n; ++i) free(walk_alns[i].target_aln);
            dk_stats_add_attempt(chunk, noisy_reg_beg, noisy_reg_end, exp_beg, exp_end, n_reads, n_noisy_reads, subsampled, n_walks, walks_null, demoted, sel, collapsed, asm_secs, pass, !junc_ok ? DK_OUT_PAIR_NO_JUNC : rescue_n == 0 ? DK_OUT_POOR_ALN : DK_OUT_PAIR_NO_SV);
            dk_walks_free(walks);
            for (int i = 0; i < n_reads; ++i) { free(read_strs[i]); free(read_names[i]); }
            free(read_strs); free(read_names); free(read_ids); free(ref_seq);
            return -1;
        }
    }
    if (aln_ok) {
        hts_pos_t win_beg = reg_beg - opt->noisy_reg_flank_len, win_end = reg_end + opt->noisy_reg_flank_len;
        if (paired) {
            // the variant window covers both breakpoints and the predicted SV with a flank long enough for the
            // read-vs-consensus alignments to reach through the SV junction (the ordinary 10 bp flank leaves a
            // consensus that ends a few bases past the junction, and reads then do not count as covering it)
            win_beg = (pred_sv_beg >= 0 && pred_sv_beg < reg_beg ? pred_sv_beg : reg_beg) - LONGCALLD_DEKNOT_PAIR_TRIM_FLANK;
            win_end = (pred_sv_end > reg_end ? pred_sv_end : reg_end) + LONGCALLD_DEKNOT_PAIR_TRIM_FLANK;
        }
        if (win_beg < exp_beg) win_beg = exp_beg;
        if (win_end > exp_end) win_end = exp_end;
        for (int iter = 0; iter < 10; ++iter) {
            int extended = 0;
            for (int i = 0; i < rescue_n; ++i) extended |= deknot_extend_ref_window(walk_alns + i, exp_beg, &win_beg, &win_end);
            if (win_beg <= exp_beg) { win_beg = exp_beg; *hit_window_edge = 1; }
            if (win_end >= exp_end) { win_end = exp_end; *hit_window_edge = 1; }
            if (!extended) break;
        }
        if (*hit_window_edge && abort_on_edge) {
            dk_stats_add_attempt(chunk, noisy_reg_beg, noisy_reg_end, exp_beg, exp_end, n_reads, n_noisy_reads, subsampled, n_walks, walks_null, demoted, sel, collapsed, asm_secs, pass, DK_OUT_HIT_EDGE);
            for (int i = 0; i < rescue_n; ++i) {
                if (walk_alns[i].target_aln != NULL) free(walk_alns[i].target_aln);
            }
            dk_walks_free(walks);
            for (int i = 0; i < n_reads; ++i) { free(read_strs[i]); free(read_names[i]); }
            free(read_strs); free(read_names); free(read_ids); free(ref_seq);
            return -1;
        }
        uint8_t *cons_seqs[2] = {NULL, NULL}; int cons_lens[2] = {0, 0};
        int trim_ok = 1;
        for (int i = 0; i < rescue_n; ++i) deknot_trim_aln(walk_alns + i, exp_beg, win_beg, win_end);
        if (paired && rescue_n == 2) {
            hts_pos_t b0 = 0, e0 = 0, b1 = 0, e1 = 0;
            int n0 = deknot_large_del_union(walk_alns + 0, win_beg, LONGCALLD_DEKNOT_PAIR_MIN_SV, &b0, &e0);
            int n1 = deknot_large_del_union(walk_alns + 1, win_beg, LONGCALLD_DEKNOT_PAIR_MIN_SV, &b1, &e1);
            if (n0 > 0 && n1 > 0 && b0 <= e1 && b1 <= e0) {
                for (int i = 0; i < 2; ++i) {
                    int merged = deknot_merge_large_dels(walk_alns + i, LONGCALLD_DEKNOT_PAIR_MIN_SV);
                    if (merged && LONGCALLD_VERBOSE >= 1)
                        fprintf(stderr, "DeKnotMergeDel %s:%" PRIi64 "-%" PRIi64 " walk=%d: large deletion pieces merged into one record\n", chunk->tname, noisy_reg_beg, noisy_reg_end, i);
                }
            }
        }
        for (int i = 0; i < rescue_n; ++i) {
            if (walk_alns[i].aln_len <= 0) trim_ok = 0;
            // trimmed walk sequence in cons space
            cons_seqs[i] = (uint8_t*)malloc((walk_alns[i].aln_len > 0 ? walk_alns[i].aln_len : 1) * sizeof(uint8_t));
            for (int c = 0; c < walk_alns[i].aln_len; ++c) {
                if (walk_alns[i].query_aln[c] != 5) cons_seqs[i][cons_lens[i]++] = walk_alns[i].query_aln[c];
            }
            if (cons_lens[i] <= 0) trim_ok = 0;
        }
        // two walks with the same sequence over the noisy region (DeKnot emits two walks per phased component
        // even when nothing separates them) would tie on every read and have one demoted: keep the
        // better-aligning one as a single consensus instead
        if (trim_ok && rescue_n == 2) {
            int gap = 0, ed = deknot_pair_diff(cons_seqs[0], cons_lens[0], cons_seqs[1], cons_lens[1], &gap);
            int max_len = (cons_lens[0] > cons_lens[1]) ? cons_lens[0] : cons_lens[1];
            double max_ed = LONGCALLD_DEKNOT_IDENT_MAX_FRAC * max_len;
            if (max_ed < LONGCALLD_DEKNOT_IDENT_MAX_ED) max_ed = LONGCALLD_DEKNOT_IDENT_MAX_ED;
            if (ed >= 0 && ed <= max_ed && gap < opt->min_sv_len) {
                int drop = (walk_score[1] > walk_score[0]) ? 0 : 1;
                free(walk_alns[drop].target_aln); free(cons_seqs[drop]);
                if (drop == 0) { walk_alns[0] = walk_alns[1]; cons_seqs[0] = cons_seqs[1]; cons_lens[0] = cons_lens[1]; walk_score[0] = walk_score[1]; }
                walk_alns[1].target_aln = NULL; walk_alns[1].query_aln = NULL; walk_alns[1].aln_len = 0;
                cons_seqs[1] = NULL; cons_lens[1] = 0;
                rescue_n = 1; collapsed = 1;
            }
        }
        if (trim_ok) { // hand the trimmed walks to the caller's slots
            for (int i = 0; i < rescue_n; ++i) {
                aln_str_t *dst = LONGCALLD_REF_CONS_ALN_STR(aln_strs[i]);
                if (dst->target_aln != NULL) free(dst->target_aln);
                *dst = walk_alns[i];
                if (clu_read_ids[i] != NULL) free(clu_read_ids[i]);
                clu_read_ids[i] = (int*)malloc(n_noisy_reads * sizeof(int));
                clu_n_seqs[i] = 0;
            }
        } else {
            for (int i = 0; i < rescue_n; ++i) free(walk_alns[i].target_aln);
        }
        // assign each noisy-region read to the walk it matches best, and record its read-vs-cons alignment
        // (clusters are restricted to the noisy-region read set: the aln_str slots are sized for it)
        double reads_t0 = realtime();
        for (int k = 0; trim_ok && k < n_noisy_reads; ++k) {
            int rlen = 0, rcover = 0, rb = 0, re = 0; uint8_t *rseq = NULL;
            collect_noisy_read_info1(opt, chunk, noisy_reads[k], win_beg, win_end, &rlen, &rseq, &rcover, &rb, &re);
            // reads fully interior to the window carry no boundary anchor
            if (rlen <= 0 || LONGCALLD_NOISY_IS_NOT_COVER(rcover)) { free(rseq); continue; }
            aln_str_t cand[2];
            int scores[2] = {INT32_MIN, INT32_MIN}, n_eqs[2] = {0, 0}, n_xids[2] = {0, 0};
            for (int i = 0; i < rescue_n; ++i) {
                cand[i].target_aln = NULL; cand[i].query_aln = NULL; cand[i].aln_len = 0;
                cand[i].target_beg = 0; cand[i].target_end = -1; cand[i].query_beg = 0; cand[i].query_end = -1;
                wfa_collect_aln_str(opt, cons_seqs[i], cons_lens[i], rseq, rlen, rcover,
                                    LONGCALLD_WFA_NO_HEURISTIC, LONGCALLD_WFA_AFFINE_2P, cand + i);
                if (cand[i].target_aln == NULL || cand[i].aln_len <= 0) continue;
                int eq = 0, xid = 0;
                // score the columns the read occupies (a partial alignment carries the rest of the consensus as gap columns)
                int c_beg = cand[i].target_beg > cand[i].query_beg ? cand[i].target_beg : cand[i].query_beg;
                int c_end = cand[i].target_end < cand[i].query_end ? cand[i].target_end : cand[i].query_end;
                for (int c = c_beg; c <= c_end && c < cand[i].aln_len; ++c) {
                    if (c < 0) continue;
                    if (cand[i].target_aln[c] == cand[i].query_aln[c] && cand[i].target_aln[c] != 5) eq++;
                    else xid++;
                }
                scores[i] = eq - xid; n_eqs[i] = eq; n_xids[i] = xid;
            }
            int best_i = (rescue_n == 2 && scores[1] > scores[0]) ? 1 : 0;
            // a read that reaches only one window edge is aligned over a different length on each walk (all of it on
            // the long one, the part that fits on a short deletion walk), so eq - xid would favor the longer alignment
            // whatever the read's haplotype: compare error rates, and prefer the longer alignment on a tie
            if (rescue_n == 2 && !LONGCALLD_NOISY_IS_BOTH_COVER(rcover) && scores[0] != INT32_MIN && scores[1] != INT32_MIN) {
                double r0 = (double)n_xids[0] / (n_eqs[0] + n_xids[0] + 1), r1 = (double)n_xids[1] / (n_eqs[1] + n_xids[1] + 1);
                if (r1 < r0) best_i = 1;
                else if (r0 < r1) best_i = 0;
                else best_i = (n_eqs[1] > n_eqs[0]) ? 1 : 0;
            }
            if (rescue_n == 2) {
                int other = 1 - best_i;
                if (cand[other].target_aln != NULL) free(cand[other].target_aln);
            }
            if (scores[best_i] == INT32_MIN || cand[best_i].target_aln == NULL) {
                if (cand[best_i].target_aln != NULL) free(cand[best_i].target_aln);
                free(rseq);
                continue;
            }
            int j = clu_n_seqs[best_i];
            aln_str_t *cr = LONGCALLD_CONS_READ_ALN_STR(aln_strs[best_i], j);
            if (cr->target_aln != NULL) free(cr->target_aln);
            *cr = cand[best_i];
            free(rseq);
            clu_read_ids[best_i][clu_n_seqs[best_i]++] = noisy_reads[k];
        }
        for (int i = 0; i < rescue_n; ++i) free(cons_seqs[i]);
        dk_read_secs = realtime() - reads_t0;
        // a committed walk with no read support is not a haplotype: demote to a single consensus
        if (trim_ok && rescue_n == 2 && (clu_n_seqs[0] == 0 || clu_n_seqs[1] == 0)) {
            if (clu_n_seqs[0] == 0) {
                aln_str_t *tmp_a = aln_strs[0]; aln_strs[0] = aln_strs[1]; aln_strs[1] = tmp_a;
                int *tmp_r = clu_read_ids[0]; clu_read_ids[0] = clu_read_ids[1]; clu_read_ids[1] = tmp_r;
                clu_n_seqs[0] = clu_n_seqs[1];
            }
            aln_str_t *drop = LONGCALLD_REF_CONS_ALN_STR(aln_strs[1]);
            if (drop->target_aln != NULL) { free(drop->target_aln); drop->target_aln = NULL; drop->query_aln = NULL; drop->aln_len = 0; }
            if (clu_read_ids[1] != NULL) { free(clu_read_ids[1]); clu_read_ids[1] = NULL; }
            clu_n_seqs[1] = 0;
            rescue_n = 1;
            demoted = 1;
        }
        if (trim_ok && rescue_n == 1 && clu_n_seqs[0] == 0) rescue_n = 0;
        if (!trim_ok) rescue_n = 0;
        // paired pass: the predicted SV has to survive trimming, collapsing and demotion in a committed walk
        int pair_no_sv = 0;
        if (paired && rescue_n >= 1) {
            hts_pos_t best_sb = -1, best_se = -1; int best_len = 0;
            for (int i = 0; i < rescue_n; ++i) {
                hts_pos_t sb = -1, se = -1;
                int len = deknot_walk_sv_span(LONGCALLD_REF_CONS_ALN_STR(aln_strs[i]), win_beg, gap_beg, gap_end, min_sv, &sb, &se);
                if (len > best_len) { best_len = len; best_sb = sb; best_se = se; }
            }
            if (best_len == 0) { rescue_n = 0; pair_no_sv = 1; }
            else if (sv_beg != NULL && sv_end != NULL) { *sv_beg = best_sb; *sv_end = best_se; }
        }
        if (pair_no_sv) outcome = DK_OUT_PAIR_NO_SV;
        else if (rescue_n == 0) outcome = DK_OUT_NO_SUPPORT;
        else if (rescue_n <= n_cons_existing) outcome = DK_OUT_NO_GAIN;
        else outcome = (rescue_n == 2) ? DK_OUT_COMMIT2 : DK_OUT_COMMIT1;
        if (rescue_n <= n_cons_existing) rescue_n = -1;
        *var_reg_beg = win_beg;
    } else {
        for (int i = 0; i < rescue_n; ++i) {
            if (walk_alns[i].target_aln != NULL) free(walk_alns[i].target_aln);
        }
        rescue_n = -1;
        outcome = DK_OUT_POOR_ALN;
    }
    dk_stats_add_attempt(chunk, noisy_reg_beg, noisy_reg_end, exp_beg, exp_end, n_reads, n_noisy_reads, subsampled, n_walks, walks_null, demoted, sel, collapsed, asm_secs, pass, outcome);
    dk_walks_free(walks);
    for (int i = 0; i < n_reads; ++i) { free(read_strs[i]); free(read_names[i]); }
    free(read_strs); free(read_names); free(read_ids); free(ref_seq);
    return rescue_n;
}

/* partner breakpoint of a noisy region: where the clipped tails of its long-clipped reads continue on the reference */
#define LONGCALLD_DEKNOT_PAIR_MAX_TAIL 5000      // at most this many tail bases (nearest the clip) are aligned
#define LONGCALLD_DEKNOT_PAIR_MAX_CLIP_READS 16  // tails placed per region
#define LONGCALLD_DEKNOT_PAIR_CLUSTER 2000       // partner positions within this distance agree
#define LONGCALLD_DEKNOT_PAIR_EXCL 2000          // the partner lies outside the noisy region +- this
#define LONGCALLD_DEKNOT_PAIR_MIN_SA_SPAN 300    // an SA entry needs this many aligned reference bases
#define LONGCALLD_DEKNOT_PAIR_MIN_READS 3        // reads that have to agree on the partner breakpoint
#define LONGCALLD_DEKNOT_PAIR_STRONG_TAIL 2000   // a placed tail at least this long counts at the read-type divergence;
#define LONGCALLD_DEKNOT_PAIR_SHORT_TAIL_DIV 0.05 // a shorter one only when it places with at most this divergence
#define LONGCALLD_DEKNOT_PAIR_UNIQUE_MARGIN 0.05 // the runner-up placement has to be worse by this fraction of the tail
#define LONGCALLD_DEKNOT_PAIR_UNIQUE_MIN_MARGIN 10 // (at least this many edits)

// first SA entry on contig tname that lies outside [reg_beg-EXCL, reg_end+EXCL] within max_dist of it:
// returns 1 and the entry's breakpoint (its end nearest the region)
static int deknot_sa_partner(const char *sa, const char *tname, hts_pos_t reg_beg, hts_pos_t reg_end, hts_pos_t max_dist, hts_pos_t *bp) {
    size_t tlen = strlen(tname);
    const char *p = sa;
    while (*p != '\0') {
        const char *end = strchr(p, ';');
        if (end == NULL) end = p + strlen(p);
        const char *c1 = memchr(p, ',', end - p);
        if (c1 != NULL && (size_t)(c1 - p) == tlen && strncmp(p, tname, tlen) == 0) {
            hts_pos_t pos = strtoll(c1 + 1, NULL, 10);
            const char *c2 = memchr(c1 + 1, ',', end - c1 - 1);          // after pos
            const char *c3 = c2 ? memchr(c2 + 1, ',', end - c2 - 1) : NULL; // after strand: CIGAR follows
            if (c2 != NULL && c3 != NULL && pos > 0) {
                hts_pos_t span = 0; long n = 0;
                for (const char *c = c3 + 1; c < end && *c != ','; ++c) {
                    if (*c >= '0' && *c <= '9') n = n * 10 + (*c - '0');
                    else { if (*c == 'M' || *c == 'D' || *c == 'N' || *c == '=' || *c == 'X') span += n; n = 0; }
                }
                if (span >= LONGCALLD_DEKNOT_PAIR_MIN_SA_SPAN) {
                    hts_pos_t sa_beg = pos, sa_end = pos + span - 1;
                    if (sa_beg > reg_end + LONGCALLD_DEKNOT_PAIR_EXCL && sa_beg - reg_end <= max_dist) { *bp = sa_beg; return 1; }
                    if (sa_end < reg_beg - LONGCALLD_DEKNOT_PAIR_EXCL && reg_beg - sa_end <= max_dist) { *bp = sa_end; return 1; }
                }
            }
        }
        p = (*end == ';') ? end + 1 : end;
    }
    return 0;
}

// place a clipped tail (upper-case bases) on chunk reference [search_beg, search_end] (1-based, inclusive) by
// semi-global edlib alignment. returns the edit distance (-1: no placement within max_ed, or an ambiguous one),
// the breakpoint (the tail's first reference base for a right clip, its last for a left clip) and in *second_ed
// the edit distance of the best placement elsewhere in the window (-1: none within max_ed + margin)
static int deknot_place_tail(bam_chunk_t *chunk, const char *tail, int tail_len, int right_clip,
                             hts_pos_t search_beg, hts_pos_t search_end, int max_ed, int margin, hts_pos_t *bp, int *second_ed) {
    *second_ed = -1;
    if (search_beg < chunk->ref_beg) search_beg = chunk->ref_beg;
    if (search_end > chunk->ref_end) search_end = chunk->ref_end;
    if (search_end - search_beg + 1 < tail_len) return -1;
    int tlen = (int)(search_end - search_beg + 1);
    char *target = (char*)malloc(tlen);
    for (int i = 0; i < tlen; ++i) target[i] = toupper(chunk->ref_seq[search_beg + i - chunk->ref_beg]);
    EdlibAlignResult r = edlibAlign(tail, tail_len, target, tlen, edlibNewAlignConfig(max_ed, EDLIB_MODE_HW, EDLIB_TASK_LOC, NULL, 0));
    int ed = -1, mask_lo = -1, mask_hi = -1;
    if (r.status == EDLIB_STATUS_OK && r.editDistance >= 0 && r.editDistance <= max_ed && r.numLocations >= 1
        && r.startLocations != NULL && r.endLocations != NULL) {
        // every equally good placement has to agree on the breakpoint
        int lo = right_clip ? r.startLocations[0] : r.endLocations[0], hi = lo;
        mask_lo = r.startLocations[0]; mask_hi = r.endLocations[0];
        for (int k = 1; k < r.numLocations; ++k) {
            int v = right_clip ? r.startLocations[k] : r.endLocations[k];
            if (v < lo) lo = v;
            if (v > hi) hi = v;
            if (r.startLocations[k] < mask_lo) mask_lo = r.startLocations[k];
            if (r.endLocations[k] > mask_hi) mask_hi = r.endLocations[k];
        }
        if (hi - lo <= LONGCALLD_DEKNOT_PAIR_CLUSTER) { ed = r.editDistance; *bp = search_beg + lo; }
    }
    edlibFreeAlignResult(r);
    if (ed >= 0) { // best placement elsewhere: mask the found one (plus the cluster distance) and search again
        mask_lo -= LONGCALLD_DEKNOT_PAIR_CLUSTER; mask_hi += LONGCALLD_DEKNOT_PAIR_CLUSTER;
        if (mask_lo < 0) mask_lo = 0;
        if (mask_hi > tlen - 1) mask_hi = tlen - 1;
        for (int i = mask_lo; i <= mask_hi; ++i) target[i] = 'N';
        EdlibAlignResult r2 = edlibAlign(tail, tail_len, target, tlen, edlibNewAlignConfig(ed + margin, EDLIB_MODE_HW, EDLIB_TASK_DISTANCE, NULL, 0));
        if (r2.status == EDLIB_STATUS_OK && r2.editDistance >= 0) *second_ed = r2.editDistance;
        edlibFreeAlignResult(r2);
    }
    free(target);
    return ed;
}

// find the partner breakpoint of a noisy region from its reads with long clips at the region: the SA tag when
// one points at the same contig within LONGCALLD_DEKNOT_PAIR_MAX_DIST (weight 2), otherwise the clipped tail placed
// on the chunk reference on either side of the region (weight 1; the placement has to be unique within the search
// distance by a clear margin, and a tail shorter than LONGCALLD_DEKNOT_PAIR_STRONG_TAIL only counts when it places
// with at most LONGCALLD_DEKNOT_PAIR_SHORT_TAIL_DIV divergence). returns the number of reads agreeing on the partner
// (0, or >= LONGCALLD_DEKNOT_PAIR_MIN_READS) and the partner interval [pair_beg, pair_end] spanned by their breakpoints.
int deknot_find_partner(const call_var_opt_t *opt, bam_chunk_t *chunk, hts_pos_t noisy_reg_beg, hts_pos_t noisy_reg_end,
                        int n_noisy_reads, int *noisy_reads, hts_pos_t *pair_beg, hts_pos_t *pair_end) {
    *pair_beg = -1; *pair_end = -1;
    if (chunk->reads == NULL) return 0;
    hts_pos_t bps[LONGCALLD_DEKNOT_PAIR_MAX_CLIP_READS]; int bp_w[LONGCALLD_DEKNOT_PAIR_MAX_CLIP_READS]; int n_bp = 0;
    hts_pos_t max_dist = LONGCALLD_DEKNOT_PAIR_MAX_DIST;
    double max_div = opt->is_ont ? 0.25 : 0.10;
    for (int k = 0; k < n_noisy_reads && n_bp < LONGCALLD_DEKNOT_PAIR_MAX_CLIP_READS; ++k) {
        int rid = noisy_reads[k];
        digar_t *d = chunk->digars + rid;
        if (d->n_digar <= 0) continue;
        digar1_t *first = d->digars, *last = d->digars + d->n_digar - 1;
        int right_clip = -1, clip_len = 0, clip_type = 0;
        if ((first->type == BAM_CSOFT_CLIP || first->type == BAM_CHARD_CLIP) && first->len >= LONGCALLD_DEKNOT_PAIR_MIN_CLIP
            && d->beg >= noisy_reg_beg - 2000 && d->beg <= noisy_reg_end + 2000) { right_clip = 0; clip_len = first->len; clip_type = first->type; }
        else if ((last->type == BAM_CSOFT_CLIP || last->type == BAM_CHARD_CLIP) && last->len >= LONGCALLD_DEKNOT_PAIR_MIN_CLIP
            && d->end >= noisy_reg_beg - 2000 && d->end <= noisy_reg_end + 2000) { right_clip = 1; clip_len = last->len; clip_type = last->type; }
        if (right_clip < 0) continue;
        bam1_t *b = chunk->reads[rid];
        if (b == NULL) continue;
        hts_pos_t bp = -1;
        char *sa = get_aux_str_from_bam(b, "SA"); // raw aux field: type byte, then the string
        if (sa != NULL && sa[0] == 'Z' && deknot_sa_partner(sa + 1, chunk->tname, noisy_reg_beg, noisy_reg_end, max_dist, &bp)) {
            if (LONGCALLD_VERBOSE >= 2)
                fprintf(stderr, "DeKnotTail %s:%" PRIi64 "-%" PRIi64 " %s %s_clip=%d SA partner=%" PRIi64 "\n", chunk->tname, noisy_reg_beg, noisy_reg_end,
                        bam_get_qname(b), right_clip ? "right" : "left", clip_len, bp);
            bp_w[n_bp] = 2; bps[n_bp++] = bp;
            continue;
        }
        // no usable SA entry (minimap2 leaves paralogous tails unmapped or as secondary alignments): place the tail ourselves
        if (clip_type != BAM_CSOFT_CLIP) continue;
        int qlen = b->core.l_qseq;
        if (clip_len > qlen) continue;
        int tail_len = clip_len < LONGCALLD_DEKNOT_PAIR_MAX_TAIL ? clip_len : LONGCALLD_DEKNOT_PAIR_MAX_TAIL;
        int q0 = right_clip ? qlen - clip_len : clip_len - tail_len; // tail bases nearest the clip
        char *tail = (char*)malloc(tail_len);
        uint8_t *bseq = bam_get_seq(b);
        for (int j = 0; j < tail_len; ++j) tail[j] = seq_nt16_str[bam_seqi(bseq, q0 + j)];
        hts_pos_t clip_pos = right_clip ? d->end : d->beg;
        int max_ed = (int)((tail_len < LONGCALLD_DEKNOT_PAIR_STRONG_TAIL ? LONGCALLD_DEKNOT_PAIR_SHORT_TAIL_DIV : max_div) * tail_len);
        int margin = (int)(LONGCALLD_DEKNOT_PAIR_UNIQUE_MARGIN * tail_len);
        if (margin < LONGCALLD_DEKNOT_PAIR_UNIQUE_MIN_MARGIN) margin = LONGCALLD_DEKNOT_PAIR_UNIQUE_MIN_MARGIN;
        // reads are only loaded for the chunk region: search within it, on either side of the noisy region
        hts_pos_t lo = chunk->reg_beg > clip_pos - max_dist ? chunk->reg_beg : clip_pos - max_dist;
        hts_pos_t hi = chunk->reg_end < clip_pos + max_dist ? chunk->reg_end : clip_pos + max_dist;
        hts_pos_t lbp = -1, rbp = -1; int lsec = -1, rsec = -1;
        int led = deknot_place_tail(chunk, tail, tail_len, right_clip, lo, noisy_reg_beg - LONGCALLD_DEKNOT_PAIR_EXCL, max_ed, margin, &lbp, &lsec);
        int red = deknot_place_tail(chunk, tail, tail_len, right_clip, noisy_reg_end + LONGCALLD_DEKNOT_PAIR_EXCL, hi, max_ed, margin, &rbp, &rsec);
        free(tail);
        // best placement over both sides; the runner-up is the other side's best or the same side's second
        int ed = -1, second = -1;
        if (led >= 0 && (red < 0 || led <= red)) { ed = led; bp = lbp; second = lsec; if (red >= 0 && (second < 0 || red < second)) second = red; }
        else if (red >= 0) { ed = red; bp = rbp; second = rsec; if (led >= 0 && (second < 0 || led < second)) second = led; }
        int unique = (ed >= 0 && (second < 0 || second > ed + margin));
        if (LONGCALLD_VERBOSE >= 2)
            fprintf(stderr, "DeKnotTail %s:%" PRIi64 "-%" PRIi64 " %s %s_clip=%d tail=%d placed=%" PRIi64 " ed=%d second=%d margin=%d unique=%d\n",
                    chunk->tname, noisy_reg_beg, noisy_reg_end, bam_get_qname(b), right_clip ? "right" : "left", clip_len, tail_len, bp, ed, second, margin, unique);
        if (unique) { bp_w[n_bp] = 1; bps[n_bp++] = bp; }
    }
    if (n_bp < LONGCALLD_DEKNOT_PAIR_MIN_READS) return 0;
    for (int i = 1; i < n_bp; ++i) { // insertion sort by position
        hts_pos_t v = bps[i]; int w = bp_w[i]; int j = i;
        while (j > 0 && bps[j-1] > v) { bps[j] = bps[j-1]; bp_w[j] = bp_w[j-1]; j--; }
        bps[j] = v; bp_w[j] = w;
    }
    // best cluster: the heaviest (SA partners count double), then the largest
    int best_i = -1, best_j = -1, best_w = 0;
    for (int i = 0; i < n_bp; ++i) {
        int j = i, w = bp_w[i];
        while (j + 1 < n_bp && bps[j+1] - bps[i] <= LONGCALLD_DEKNOT_PAIR_CLUSTER) { j++; w += bp_w[j]; }
        if (j - i + 1 < LONGCALLD_DEKNOT_PAIR_MIN_READS) continue;
        if (w > best_w || (w == best_w && j - i > best_j - best_i)) { best_i = i; best_j = j; best_w = w; }
    }
    if (best_i < 0) return 0;
    *pair_beg = bps[best_i]; *pair_end = bps[best_j];
    dk_stats_count(&dk_stats.n_pair_regions);
    return best_j - best_i + 1;
}

// local-assembly rescue over one window covering the noisy region and its partner breakpoint [pair_beg, pair_end]
// (reads: both breakpoint regions' reads). returns the committed consensus count, or -1; on a commit the walks carry
int deknot_sv_rescue_paired(const call_var_opt_t *opt, bam_chunk_t *chunk, hts_pos_t noisy_reg_beg, hts_pos_t noisy_reg_end,
                            hts_pos_t pair_beg, hts_pos_t pair_end, int n_noisy_reads, int *noisy_reads,
                            int n_cons_existing, int *clu_n_seqs, int **clu_read_ids, aln_str_t **aln_strs, hts_pos_t *var_reg_beg,
                            hts_pos_t *sv_beg, hts_pos_t *sv_end, int depth_triggered, hts_pos_t win_lo, hts_pos_t win_hi, int *junc_reads) {
    int hit_edge = 0;
    *sv_beg = -1; *sv_end = -1;
    int ret = deknot_sv_rescue_core(opt, chunk, noisy_reg_beg, noisy_reg_end, pair_beg, pair_end, n_noisy_reads, noisy_reads,
                                    LONGCALLD_DEKNOT_PAIR_FLANK_LEN, 0, &hit_edge, depth_triggered ? DK_PASS_DEPTH : DK_PASS_PAIRED,
                                    n_cons_existing, clu_n_seqs, clu_read_ids, aln_strs, var_reg_beg, sv_beg, sv_end, win_lo, win_hi, junc_reads);
    dk_stats_add_region(n_cons_existing, ret);
    return ret;
}

void deknot_stats_add_depth_runs(int n) {
    pthread_mutex_lock(&dk_stats.lock);
    dk_stats.n_depth_runs += n;
    pthread_mutex_unlock(&dk_stats.lock);
}

int deknot_sv_rescue(const call_var_opt_t *opt, bam_chunk_t *chunk, hts_pos_t noisy_reg_beg, hts_pos_t noisy_reg_end,
                     int n_noisy_reads, int *noisy_reads,
                     int n_cons_existing, int *clu_n_seqs, int **clu_read_ids, aln_str_t **aln_strs, hts_pos_t *var_reg_beg) {
    int hit_edge = 0;
    // long clips at the region are the signature of an SV larger than the small window
    int n_long_clip = 0;
    for (int k = 0; k < n_noisy_reads; ++k) {
        digar_t *d = chunk->digars + noisy_reads[k];
        if (d->n_digar <= 0) continue;
        digar1_t *first = d->digars, *last = d->digars + d->n_digar - 1;
        if ((first->type == BAM_CSOFT_CLIP || first->type == BAM_CHARD_CLIP) && first->len >= 500
            && d->beg >= noisy_reg_beg - 2000 && d->beg <= noisy_reg_end + 2000) n_long_clip++;
        else if ((last->type == BAM_CSOFT_CLIP || last->type == BAM_CHARD_CLIP) && last->len >= 500
            && d->end >= noisy_reg_beg - 2000 && d->end <= noisy_reg_end + 2000) n_long_clip++;
    }
    int ret;
    if (n_long_clip >= 2) {
        ret = deknot_sv_rescue_core(opt, chunk, noisy_reg_beg, noisy_reg_end, -1, -1, n_noisy_reads, noisy_reads,
                                    LONGCALLD_DEKNOT_FLANK_LEN, 0, &hit_edge, DK_PASS_WIDE,
                                    n_cons_existing, clu_n_seqs, clu_read_ids, aln_strs, var_reg_beg, NULL, NULL, -1, -1, NULL);
        dk_stats_add_region(n_cons_existing, ret);
        return ret;
    }
    ret = deknot_sv_rescue_core(opt, chunk, noisy_reg_beg, noisy_reg_end, -1, -1, n_noisy_reads, noisy_reads,
                                LONGCALLD_DEKNOT_SMALL_FLANK_LEN, 1, &hit_edge, DK_PASS_SMALL,
                                n_cons_existing, clu_n_seqs, clu_read_ids, aln_strs, var_reg_beg, NULL, NULL, -1, -1, NULL);
    // an indel run crossing the small window's edge means an SV larger than the flank: retry wide
    if (ret < 0 && hit_edge) {
        ret = deknot_sv_rescue_core(opt, chunk, noisy_reg_beg, noisy_reg_end, -1, -1, n_noisy_reads, noisy_reads,
                                    LONGCALLD_DEKNOT_FLANK_LEN, 0, &hit_edge, DK_PASS_WIDE_RETRY,
                                    n_cons_existing, clu_n_seqs, clu_read_ids, aln_strs, var_reg_beg, NULL, NULL, -1, -1, NULL);
    }
    dk_stats_add_region(n_cons_existing, ret);
    return ret;
}

int collect_noisy_reg_aln_strs(const call_var_opt_t *opt, bam_chunk_t *chunk, hts_pos_t noisy_reg_beg, hts_pos_t noisy_reg_end, int noisy_reg_i,
                               int n_noisy_reg_reads, int *noisy_read_ids, uint8_t *ref_seq, int ref_seq_len,
                               int *clu_n_seqs, int **clu_read_ids, aln_str_t **aln_strs) {
    if (n_noisy_reg_reads <= 0) return 0;
    // fully_cover: 0 -> none, 1 -> left, 2 -> right, 3 -> both
    char **names = NULL; uint8_t **seqs = NULL; uint8_t *strands=NULL; int *fully_covers = NULL, *lens = NULL, *haps = NULL; hts_pos_t *phase_sets = NULL; uint8_t **base_quals = NULL;
    int *read_id_to_full_covers = NULL, *read_reg_beg = NULL, *read_reg_end = NULL;
    collect_noisy_read_info(opt, chunk, noisy_reg_beg, noisy_reg_end, noisy_reg_i, n_noisy_reg_reads, noisy_read_ids,
                            &lens, &seqs, &strands, &base_quals, &names, &fully_covers, &read_id_to_full_covers, 
                            &read_reg_beg, &read_reg_end, &haps, &phase_sets);

    // for large noisy regions, sort reads by fully_cover, error_base_rate, and length
    int sampling_reads = 0;
    if (noisy_reg_end - noisy_reg_beg + 1 >= opt->min_noisy_reg_size_to_sample_reads) sampling_reads = 1;
    sort_noisy_region_reads(n_noisy_reg_reads, noisy_read_ids, lens, seqs, base_quals, strands, fully_covers, names, haps, phase_sets, sampling_reads);
    // >= min_hap_full_read_count reads for each hap && >= min_hap_read_count reads (including not full-cover, but >= full-cover length) for each hap
    int min_hap_full_read_count = opt->min_hap_full_reads, min_hap_read_count = opt->min_hap_reads;
    int min_no_hap_full_read_count = opt->min_dp; // opt->min_no_hap_full_reads;
    // XXX find PS that fully cover the noisy region, allow 1 or 2 haps to have NO fully-covered reads, this can be resolved by POA-extension and consensus mering
    // if no valid PS was found, do haplotype-unaware local-assembly
    hts_pos_t ps_with_both_haps = collect_phase_set_with_both_haps(n_noisy_reg_reads, haps, lens, phase_sets, fully_covers, min_hap_full_read_count, min_hap_read_count);
    int n_full_reads = 0;
    for (int i = 0; i < n_noisy_reg_reads; ++i) if (LONGCALLD_NOISY_IS_BOTH_COVER(fully_covers[i])) n_full_reads++;
    int n_cons = 0;

    int collect_ref_read_aln_str = 0;
    if ((opt->refine_bam && opt->out_aln_fp != NULL) || opt->out_somatic) collect_ref_read_aln_str = 1;
    // XXX only use fully-covered reads, including cliping reads (after re-align to backbone read)
    // two cases to call consensus sequences
    if (ps_with_both_haps > 0) { // call consensus sequences for each haplotype
        n_cons = wfa_collect_noisy_aln_str_with_ps_hap(opt, sampling_reads, n_noisy_reg_reads, noisy_read_ids, lens, seqs, strands, base_quals, names, haps, phase_sets, fully_covers, ps_with_both_haps, min_hap_full_read_count, min_hap_read_count,
                                                       ref_seq, ref_seq_len, clu_n_seqs, clu_read_ids, aln_strs, collect_ref_read_aln_str);
        if (LONGCALLD_VERBOSE >= 1) fprintf(stderr, "Hap %s:%" PRIi64 "-%" PRIi64 " %" PRIi64 " %d reads (%d full) n_cons: %d\n", chunk->tname, noisy_reg_beg, noisy_reg_end, noisy_reg_end-noisy_reg_beg+1, n_noisy_reg_reads, n_full_reads, n_cons);
    // >= min_no_hap_full_read_count full reads in total
    } else if (ps_with_both_haps <= 0 && n_full_reads >= min_no_hap_full_read_count) {
        // XXX do NOT de novo abPOA for homopolymer regions
        n_cons = wfa_collect_noisy_aln_str_no_ps_hap(opt, n_noisy_reg_reads, noisy_read_ids, lens, seqs, names, fully_covers,
                                                     ref_seq, ref_seq_len, clu_n_seqs, clu_read_ids, aln_strs, collect_ref_read_aln_str);
        if (LONGCALLD_VERBOSE >= 1) fprintf(stderr, "NoHap %s:%" PRIi64 "-%" PRIi64 " %" PRIi64 " %d reads (%d full) n_cons: %d\n", chunk->tname, noisy_reg_beg, noisy_reg_end, noisy_reg_end-noisy_reg_beg+1, n_noisy_reg_reads, n_full_reads, n_cons);
    } else {
        if (LONGCALLD_VERBOSE >= 1) fprintf(stderr, "Skipped %s:%" PRIi64 "-%" PRIi64 " %" PRIi64 " %d reads (%d full)\n", chunk->tname, noisy_reg_beg, noisy_reg_end, noisy_reg_end-noisy_reg_beg+1, n_noisy_reg_reads, n_full_reads);
    }
    // update digar based on ref vs read in MSA
    if (n_cons > 0 && ((opt->refine_bam && opt->out_aln_fp != NULL) || opt->out_somatic)) {
        update_digars_from_aln_str(chunk, noisy_reg_beg, noisy_reg_end, read_id_to_full_covers, read_reg_beg, read_reg_end,
                                   aln_strs, n_cons, clu_n_seqs, clu_read_ids);
    }
    for (int i = 0; i < n_noisy_reg_reads; ++i) {
        free(seqs[i]); free(base_quals[i]);
    }
    free(names); free(strands); free(seqs); free(base_quals); free(lens); free(read_reg_beg); free(read_reg_end);
    free(fully_covers); free(read_id_to_full_covers); free(haps); free(phase_sets);
    return n_cons;
}
