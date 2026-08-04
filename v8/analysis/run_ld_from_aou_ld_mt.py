#!/usr/bin/env python3
"""
Run LD prune / LD matrix / LD score steps from a pre-generated AoU LD MatrixTable.

This script expects an LD-ready MT produced by `generate_ld_data_from_aou_vds.py`
with entry field `GT` and row fields:
    - gen_anc_ac
    - gen_anc_an
    - gen_anc_af
    - gen_anc_callrate
    - gen_anc_freq (from hl.agg.call_stats over GT)
"""

from __future__ import annotations

import argparse
import logging
import sys
from typing import Optional

import hail as hl
from hail.linalg import BlockMatrix
from hail.utils.misc import new_temp_file


logger = logging.getLogger("aou_ld_from_ld_mt")
logger.setLevel(logging.INFO)

# Hail scratch (short retention is fine; these are throwaway intermediates).
TMP_BUCKET = "gs://aou_tmp/ld_marten"
# Pipeline outputs. MUST match MY_BUCKET in generate_ld_data_from_aou_vds.py, which writes the
# LD MT this script reads by default.
MY_BUCKET = "gs://aou_tmp_30_days/ld_marten"

RARE_FREQ = 0.0005
COMMON_FREQ = 0.005

SINGLETON_AC = 1
DOUBLETON_AC = 2


def configure_logging(level: int = logging.INFO) -> None:
    """Ensure INFO logs are emitted in Dataproc driver stdout/stderr."""
    logging.basicConfig(
        level=level,
        format="%(asctime)s %(levelname)s %(name)s - %(message)s",
        stream=sys.stdout,
        force=True,
    )
    logger.setLevel(level)


def get_sample_count(mt: hl.MatrixTable, gen_anc: str) -> int:
    """Single-ancestry run: return cohort sample count for LD score adjustment."""
    n = mt.count_cols()
    logger.info("Sample count for %s: %s", gen_anc, n)
    return n


def _paths_ld_pruned(
    prefix: str,
    pop: str,
    ld_contig: str,
    r2: str,
    freq: Optional[float],
    ac_cutoff: Optional[int],
    adj: bool,
    custom_suffix: str,
) -> str:
    adj_s = "adj" if adj else "noadj"
    cs = f"_{custom_suffix}" if custom_suffix else ""
    return (
        f"{prefix}_pruned_{pop}_{ld_contig}_r2_{r2}_"
        f"{_filter_mode_token(freq, ac_cutoff)}_{adj_s}{cs}.ht"
    )


def _paths_ld_index(
    prefix: str,
    pop: str,
    ld_contig: str,
    freq: Optional[float],
    ac_cutoff: Optional[int],
    adj: bool,
    skip_ld_prune: bool,
    custom_suffix: str,
) -> str:
    adj_s = "adj" if adj else "noadj"
    prune_s = "noprune" if skip_ld_prune else "prunedrop"
    cs = f"_{custom_suffix}" if custom_suffix else ""
    return f"{prefix}_index_{pop}_{ld_contig}_{_filter_mode_token(freq, ac_cutoff)}_{adj_s}_{prune_s}{cs}.ht"


def _paths_ld_matrix(
    prefix: str,
    pop: str,
    ld_contig: str,
    freq: Optional[float],
    ac_cutoff: Optional[int],
    adj: bool,
    skip_ld_prune: bool,
    custom_suffix: str,
) -> str:
    adj_s = "adj" if adj else "noadj"
    prune_s = "noprune" if skip_ld_prune else "prunedrop"
    cs = f"_{custom_suffix}" if custom_suffix else ""
    return f"{prefix}_bm_{pop}_{ld_contig}_{_filter_mode_token(freq, ac_cutoff)}_{adj_s}_{prune_s}{cs}.bm"


def _paths_ld_scores(
    prefix: str,
    pop: str,
    ld_contig: str,
    freq: Optional[float],
    ac_cutoff: Optional[int],
    adj: bool,
    skip_ld_prune: bool,
    custom_suffix: str,
    call_rate_cutoff: float,
) -> str:
    adj_s = "adj" if adj else "noadj"
    prune_s = "noprune" if skip_ld_prune else "prunedrop"
    cs = f"_{custom_suffix}" if custom_suffix else ""
    return (
        f"{prefix}_ldscores_{pop}_{ld_contig}_{_filter_mode_token(freq, ac_cutoff)}"
        f"_{adj_s}_{prune_s}_cr{call_rate_cutoff}{cs}.ht"
    )


def _filter_mode_token(freq: Optional[float], ac_cutoff: Optional[int]) -> str:
    if freq is not None:
        return f"f{freq}"
    if ac_cutoff is not None:
        return f"ac{ac_cutoff}"
    return "nofilter"


def _paths_ld_filtered_mt(
    prefix: str,
    pop: str,
    ld_contig: str,
    freq: Optional[float],
    ac_cutoff: Optional[int],
    adj: bool,
    custom_suffix: str,
) -> str:
    adj_s = "adj" if adj else "noadj"
    cs = f"_{custom_suffix}" if custom_suffix else ""
    return f"{prefix}_filteredld_{pop}_{ld_contig}_{_filter_mode_token(freq, ac_cutoff)}_{adj_s}{cs}.mt"


def _paths_ld_matrix_input_mt(
    prefix: str,
    pop: str,
    ld_contig: str,
    freq: Optional[float],
    ac_cutoff: Optional[int],
    r2: str,
    adj: bool,
    skip_ld_prune: bool,
    custom_suffix: str,
) -> str:
    adj_s = "adj" if adj else "noadj"
    prune_s = "noprune" if skip_ld_prune else "prunedrop"
    cs = f"_{custom_suffix}" if custom_suffix else ""
    return (
        f"{prefix}_matrix_input_{pop}_{ld_contig}_r2_{r2}_"
        f"{_filter_mode_token(freq, ac_cutoff)}_{adj_s}_{prune_s}{cs}.mt"
    )


def _filter_rows_for_ld(
    mt: hl.MatrixTable,
    freq: Optional[float],
    ac_cutoff: Optional[int],
) -> hl.MatrixTable:
    """
    Filter rows to LD-eligible variants using post-ADJ cohort callstats (``gen_anc_freq``).

    Symmetric MAF band ``[freq, 1 - freq]`` matches the historical gnomAD-style LD *scores*
    step in ``OLD_generate_ld_data_from_vds.py`` (previously applied again when reading the
    index table). We apply it once here so prune, matrix, and scores share the same variant set.

    We intentionally use row-level ``gen_anc_freq`` computed from surviving cohort GT calls so
    variants with no true observed calls can be excluded before LD.
    """
    logger.info(
        "Starting _filter_rows_for_ld(freq=%s, ac_cutoff=%s)",
        freq,
        ac_cutoff,
    )
    # Behavior when getting all keys and checking was very odd? 
    # required = {"gen_anc_freq"}
    # row_keys = set(mt.row.dtype.keys())
    # if not required.issubset(row_keys):
    #     missing = sorted(required - row_keys)
    #     raise ValueError(
    #         f"MatrixTable is missing required row fields for LD filtering: {missing}"
    #     )

    af = mt.gen_anc_freq.AF[1]
    ac = mt.gen_anc_freq.AC[1]
    an = mt.gen_anc_freq.AN
    hom_alt = mt.gen_anc_freq.homozygote_count[1]
    cond_expr = (an - ac > 1) & ~((af == 0.5) & (hom_alt == 0))
    if ac_cutoff is not None:
        cond_expr = cond_expr & (ac > ac_cutoff)
    if freq is not None:
        cond_expr = cond_expr & (af >= freq) & (af <= 1 - freq)
    return mt.filter_rows(cond_expr)


def filter_mt_for_ld(
    mt: hl.MatrixTable,
    freq: Optional[float] = None,
    ld_pruned_path: Optional[str] = None,
    ac_cutoff: Optional[int] = None,
) -> hl.MatrixTable:
    logger.info("Starting filter_mt_for_ld")
    logger.info(
        "filter_mt_for_ld: freq=%s ac_cutoff=%s ld_pruned_path=%s",
        freq,
        ac_cutoff,
        ld_pruned_path,
    )

    if ld_pruned_path:
        # ``hl.ld_prune`` returns a maximal LD-*independent* set (one representative per correlated
        # cluster). We **exclude** those representatives from the LD matrix so downstream LD uses
        # variants that remain in dense / correlated regions (related variants), not the tag SNPs.
        logger.info("Excluding LD-prune independent representatives from %s", ld_pruned_path)
        ld_ht = hl.read_table(ld_pruned_path)
        mt = mt.filter_rows(~hl.is_defined(ld_ht[mt.row_key]))

    mt = _filter_rows_for_ld(
        mt,
        freq=freq,
        ac_cutoff=ac_cutoff,
    )
    logger.info("filter_mt_for_ld: %s variants, %s samples", *mt.count())
    return mt


def generate_ld_pruned_set(
    mt: hl.MatrixTable,
    gen_anc: str,
    r2: str,
    freq: Optional[float],
    radius: int,
    overwrite: bool,
    ld_contig: str,
    ac_cutoff: Optional[int],
    adj: bool,
    custom_suffix: str,
    output_prefix: str,
) -> None:
    logger.info(
        "Starting generate_ld_pruned_set(gen_anc=%s, ld_contig=%s, r2=%s, freq=%s, radius=%s)",
        gen_anc,
        ld_contig,
        r2,
        freq,
        radius,
    )
    # ``mt`` is expected to be pre-filtered once per contig in ``main``.
    logger.info("Select rows, cols, and globals, then massively naive coelasce to 25 partitions, before ld_prune")
    mt_temp = mt.select_rows().select_cols().select_globals()
    mt_temp = mt_temp.naive_coalesce(25).checkpoint(new_temp_file("ld_set_super_coalesce", "mt"))
    pruned_ht = hl.ld_prune(mt_temp.GT, r2=float(r2), bp_window_size=radius)
    ht = mt_temp.rows()
    # Rows in ``pruned_ht`` are LD-*independent* representatives; we persist them so the matrix
    # step can **drop** those rows and retain correlated (related) variants for LD.
    ht = ht.filter(hl.is_defined(pruned_ht[ht.key]))
    out = _paths_ld_pruned(
        output_prefix, gen_anc, ld_contig, r2, freq, ac_cutoff, adj, custom_suffix
    )
    logger.info("Writing LD pruned set to %s", out)
    ht.write(out, overwrite)


def generate_ld_matrix(
    mt: hl.MatrixTable,
    gen_anc: str,
    radius: int,
    freq: Optional[float],
    adj: bool,
    overwrite: bool,
    ld_contig: str,
    r2: str,
    custom_suffix: str,
    ac_cutoff: Optional[int],
    skip_ld_prune: bool,
    output_prefix: str,
    ld_block_size: int = 2048,
) -> None:
    logger.info(
        "Starting generate_ld_matrix(gen_anc=%s, ld_contig=%s, r2=%s, freq=%s, radius=%s, "
        "ld_block_size=%s)",
        gen_anc,
        ld_contig,
        r2,
        freq,
        radius,
        ld_block_size,
    )
    # ``mt`` is expected to be pre-filtered once per contig in ``main``.
    if skip_ld_prune:
        logger.info("Skipping LD-prune row exclusion before matrix construction (--skip-ld-prune)")
    else:
        pruned_path = _paths_ld_pruned(
            output_prefix, gen_anc, ld_contig, r2, freq, ac_cutoff, adj, custom_suffix
        )
        ld_ht = hl.read_table(pruned_path)
        mt = mt.filter_rows(~hl.is_defined(ld_ht[mt.row_key]))
    logger.info("Checkpoint MT before ld_matrix")
    mt = mt.checkpoint(
        _paths_ld_matrix_input_mt(
            output_prefix, gen_anc, ld_contig, freq, ac_cutoff, r2, adj, skip_ld_prune, custom_suffix
        ),
        overwrite=overwrite,
    )
    idx_path = _paths_ld_index(
        output_prefix, gen_anc, ld_contig, freq, ac_cutoff, adj, skip_ld_prune, custom_suffix
    )
    logger.info("Writing LD variant index %s", idx_path)
    mt.rows().add_index().write(idx_path, overwrite)

    ld = hl.ld_matrix(
        mt.GT.n_alt_alleles(), mt.locus, radius, block_size=ld_block_size
    )
    ld = ld.sparsify_triangle()
    bm_path = _paths_ld_matrix(
        output_prefix, gen_anc, ld_contig, freq, ac_cutoff, adj, skip_ld_prune, custom_suffix
    )
    logger.info("Writing LD BlockMatrix %s", bm_path)
    ld.write(bm_path, overwrite)


def compute_and_annotate_ld_score(
    ht: hl.Table, r2_adj: BlockMatrix, radius: int, out_name: str, overwrite: bool
) -> None:
    logger.info(
        "Starting compute_and_annotate_ld_score(out_name=%s, radius=%s)",
        out_name,
        radius,
    )
    starts_and_stops = hl.linalg.utils.locus_windows(ht.locus, radius, _localize=False)
    r2_adj = r2_adj._sparsify_row_intervals_expr(starts_and_stops, blocks_only=False)
    l2row = r2_adj.sum(axis=0).T
    l2col = r2_adj.sum(axis=1)
    l2 = l2row + l2col - 1
    l2_bm_tmp = new_temp_file()
    l2_tsv_tmp = new_temp_file()
    l2.write(l2_bm_tmp, force_row_major=True)
    BlockMatrix.export(l2_bm_tmp, l2_tsv_tmp)
    ht_scores = hl.import_table(l2_tsv_tmp, no_header=True, impute=True)
    ht_scores = ht_scores.add_index().rename({"f0": "ld_score"})
    ht_scores = ht_scores.key_by("idx")
    ht = ht.annotate(**ht_scores[ht.new_idx]).select_globals()
    ht.filter(hl.is_defined(ht.ld_score)).write(out_name, overwrite)


def generate_ld_scores_from_ld_matrix(
    gen_anc: str,
    n_samples: int,
    freq: Optional[float],
    ac_cutoff: Optional[int],
    call_rate_cutoff: float,
    adj: bool,
    skip_ld_prune: bool,
    radius: int,
    overwrite: bool,
    ld_contig: str,
    custom_suffix: str,
    output_prefix: str,
) -> None:
    logger.info(
        "Starting generate_ld_scores_from_ld_matrix(gen_anc=%s, ld_contig=%s, n_samples=%s, "
        "freq=%s, radius=%s)",
        gen_anc,
        ld_contig,
        n_samples,
        freq,
        radius,
    )
    if n_samples < 3:
        raise ValueError(
            f"LD score adjustment uses (n-1)/(n-2) with n=columns; need n>=3, got n={n_samples}. "
            "Use a larger cohort or run with --generate-ld-pruned-set / --generate-ld-matrix only."
        )

    idx_path = _paths_ld_index(
        output_prefix, gen_anc, ld_contig, freq, ac_cutoff, adj, skip_ld_prune, custom_suffix
    )
    ht = hl.read_table(idx_path)

    row_keys = set(ht.row.dtype.keys())
    if "idx" not in row_keys:
        raise ValueError("LD index HT must contain an 'idx' column from add_index().")

    # ``idx`` = matrix row id when the LD BlockMatrix was built (use for ``BlockMatrix.filter``).
    # ``new_idx`` = dense 0..n-1 over the current table after any optional future filters; use in
    # ``compute_and_annotate_ld_score`` to align exported score rows. Apply any future ``ht.filter``
    # *before* ``add_index(new_idx)``. No extra filters here yet.
    ht = ht.add_index(name="new_idx")
    indices = ht.idx.collect()

    bm_path = _paths_ld_matrix(
        output_prefix, gen_anc, ld_contig, freq, ac_cutoff, adj, skip_ld_prune, custom_suffix
    )
    r2 = BlockMatrix.read(bm_path)
    r2 = r2.filter(indices, indices) ** 2
    r2_adj = ((n_samples - 1.0) / (n_samples - 2.0)) * r2 - (1.0 / (n_samples - 2.0))

    out_name = _paths_ld_scores(
        output_prefix, gen_anc, ld_contig, freq, ac_cutoff, adj, skip_ld_prune, custom_suffix, call_rate_cutoff
    )
    logger.info("Writing LD scores to %s", out_name)
    compute_and_annotate_ld_score(ht, r2_adj, radius, out_name, overwrite)


def main(args) -> None:
    configure_logging()
    gen_anc = args.gen_anc
    contig = args.contig
    if args.test and not contig:
        contig = "chr22"

    if args.batch:
        hl.init(
            backend="batch",
            tmp_dir=args.tmp_dir,
            driver_memory="highmem",
            driver_cores=8,
            worker_memory="standard",
            worker_cores=8,
            default_reference="GRCh38",
            log="/aou_ld_from_ld_mt.log",
            app_name="aou_ld_from_ld_mt",
            gcs_requester_pays_configuration=args.requester_pays,
        )
    elif args.dataproc:
        hl.init(
            tmp_dir=args.tmp_dir,
            default_reference="GRCh38",
            app_name="aou_ld_from_ld_mt_sas22",
            gcs_requester_pays_configuration=args.requester_pays,
            backend="spark",
        )
    else:
        raise ValueError("Must run in batch or dataproc mode. Please set --batch or --dataproc.")

    hl._set_flags(use_new_shuffle="1")

    output_prefix = args.output_prefix or f"{MY_BUCKET}/{gen_anc}_{contig or 'genome'}_vdsld"
    custom_suffix = args.custom_suffix or ""
    suffix_token = f"_{custom_suffix}" if custom_suffix else ""
    ld_mt_path = args.ld_mt_path or f"{output_prefix}_ld_cohort{suffix_token}.mt"
    if args.coalesce is not None and args.repartition is not None:
        raise ValueError("Please set only one of --coalesce or --repartition.")
    if args.ld_freq is None and args.ld_ac is None:
        raise ValueError("Set at least one of --ld-freq or --ld-ac.")

    logger.info("Reading LD MatrixTable from %s", ld_mt_path)
    mt = hl.read_matrix_table(ld_mt_path)

    run_ld_steps = (
        args.generate_ld_pruned_set
        or args.generate_ld_matrix
        or args.generate_ld_scores
    )
    if not run_ld_steps:
        logger.info("LD prune/matrix/scores not requested; done.")
        return

    n_samples = get_sample_count(mt, gen_anc)

    # Fixed chromosome order only (no ``aggregate_rows`` / ``collect_as_set`` over the full MT:
    # that is too slow at AoU scale). Pass ``--contig`` when the MT is single-chromosome or you
    # want to avoid iterating absent chromosomes.
    contig_list = (
        [f"chr{i}" for i in range(1, 23)] + ["chrX", "chrY"] if not contig else [contig]
    )

    for ld_contig in contig_list:
        mt_ld = mt.filter_rows(mt.locus.contig == ld_contig)

        if args.coalesce is not None:
            logger.info(
                "Naively coalescing LD MT for %s to %s partitions",
                ld_contig,
                args.coalesce,
            )
            mt_ld = mt_ld.naive_coalesce(args.coalesce)
            mt_ld = mt_ld.checkpoint(new_temp_file("naive_coalesced", "mt"))
        elif args.repartition is not None:
            logger.info(
                "Repartitioning LD MT for %s to %s partitions",
                ld_contig,
                args.repartition,
            )
            # Shuffle is needed to enforce uniform partitions. 
            mt_ld = mt_ld.repartition(args.repartition, shuffle=True)
            mt_ld = mt_ld.checkpoint(new_temp_file("repartitioned", "mt"))


        if mt_ld.count_rows() == 0:
            logger.info("Skipping %s: no variants in LD MT for this contig.", ld_contig)
            continue

        # Bool for whether we need to generate the filtered LD MT
        need_ld_filtered_mt = args.generate_ld_pruned_set or args.generate_ld_matrix
        mt_ld_filtered = None
        if need_ld_filtered_mt:
            logger.info("Filtering LD MT once for prune/matrix reuse on %s", ld_contig)
            # Perform the act of filtering the LD MT for the first time.
            mt_ld_filtered = filter_mt_for_ld(
                mt=mt_ld,
                freq=args.ld_freq,
                ld_pruned_path=None,
                ac_cutoff=args.ld_ac,
            )
            # Path for the filtered LD MT
            filtered_mt_path = _paths_ld_filtered_mt(
                output_prefix, gen_anc, ld_contig, args.ld_freq, args.ld_ac, args.adj, custom_suffix
            )
            logger.info("Checkpoint filtered LD MT to %s", filtered_mt_path)
            mt_ld_filtered = mt_ld_filtered.checkpoint(filtered_mt_path, overwrite=args.overwrite)

        if args.skip_ld_prune and args.generate_ld_pruned_set:
            logger.info(
                "Skipping --generate-ld-pruned-set for %s due to --skip-ld-prune",
                ld_contig,
            )
        elif args.generate_ld_pruned_set:
            generate_ld_pruned_set(
                mt=mt_ld_filtered,
                gen_anc=gen_anc,
                r2=args.r2,
                freq=args.ld_freq,
                radius=args.radius,
                overwrite=args.overwrite,
                ld_contig=ld_contig,
                ac_cutoff=args.ld_ac,
                adj=args.adj,
                custom_suffix=custom_suffix,
                output_prefix=output_prefix,
            )
        if args.generate_ld_matrix:
            generate_ld_matrix(
                mt=mt_ld_filtered,
                gen_anc=gen_anc,
                radius=args.radius,
                freq=args.ld_freq,
                adj=args.adj,
                overwrite=args.overwrite,
                ld_contig=ld_contig,
                r2=args.r2,
                custom_suffix=custom_suffix,
                ac_cutoff=args.ld_ac,
                skip_ld_prune=args.skip_ld_prune,
                output_prefix=output_prefix,
                ld_block_size=args.ld_block_size,
            )
        if args.generate_ld_scores:
            generate_ld_scores_from_ld_matrix(
                gen_anc=gen_anc,
                n_samples=n_samples,
                freq=args.ld_freq,
                ac_cutoff=args.ld_ac,
                call_rate_cutoff=args.min_call_rate,
                adj=args.adj,
                skip_ld_prune=args.skip_ld_prune,
                radius=args.radius,
                overwrite=args.overwrite,
                ld_contig=ld_contig,
                custom_suffix=custom_suffix,
                output_prefix=output_prefix,
            )


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="AoU LD: run prune/matrix/scores from a precomputed LD MT")
    parser.add_argument("--gen-anc", default="eas", help="Cohort ancestry label for naming and score correction")
    parser.add_argument("--contig", default=None, help="Restrict to one chrom (e.g. chr22)")
    parser.add_argument(
        "--output-prefix",
        default=None,
        help="Prefix for pruned/index/bm/ldscores paths and default LD MT name",
    )
    parser.add_argument(
        "--ld-mt-path",
        default=None,
        help="Path to pre-generated cohort LD MatrixTable",
    )
    parser.add_argument(
        "--min-call-rate",
        dest="min_call_rate",
        type=float,
        default=0.90,
        help=(
            "Used in LD-score output path tokens only; row filtering uses symmetric MAF and AC."
        ),
    )
    parser.add_argument(
        "--generate-ld-pruned-set",
        action="store_true",
        help="Write LD-pruned variant HT (before --generate-ld-matrix)",
    )
    parser.add_argument(
        "--skip-ld-prune",
        action="store_true",
        help="Skip hl.ld_prune and do not require/read LD-pruned HT for matrix construction.",
    )
    parser.add_argument("--generate-ld-matrix", action="store_true")
    parser.add_argument("--generate-ld-scores", action="store_true")
    parser.add_argument("--r2", default="0.2", help="ld_prune r^2 threshold")
    parser.add_argument("--radius", type=int, default=1_000_000, help="LD window (bp)")
    parser.add_argument(
        "--ld-freq",
        type=float,
        default=None,
        help="AF lower bound for LD filters (can be combined with --ld-ac).",
    )
    parser.add_argument(
        "--ld-ac",
        type=int,
        default=None,
        help="AC lower bound for LD filters (can be combined with --ld-freq).",
    )
    
    parser.add_argument(
        "--adj",
        action="store_true",
        help="Tag outputs as adj-mode (input MT should already reflect adj filtering).",
    )
    parser.add_argument(
        "--ld-block-size",
        type=int,
        default=2048,
        help="hl.ld_matrix block_size",
    )
    parser.add_argument("--custom-suffix", default="", help="Suffix token in output paths")
    parser.add_argument("--overwrite", action="store_true")
    parser.add_argument("--dataproc", action="store_true")
    parser.add_argument("--batch", action="store_true", help="Run in QueryOnBatch, not Dataproc.")
    parser.add_argument("--test", action="store_true", help="Use chr22 if --contig not set")
    parser.add_argument("--tmp-dir", default=TMP_BUCKET)
    parser.add_argument(
        "--requester-pays",
        default="aou-neale-gwas",
        help="GCS requester pays project",
    )
    parser.add_argument(
        "--coalesce",
        type=int,
        default=None,
        help="Naively coalesce per-contig LD MT to N partitions before downstream steps.",
    )
    parser.add_argument(
        "--repartition",
        type=int,
        default=None,
        help="Repartition per-contig LD MT to N partitions before downstream steps.",
    )
    cli_args = parser.parse_args()
    main(cli_args)
