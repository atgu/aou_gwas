#!/usr/bin/env python3
"""
All of Us LD pipeline from a Hail VDS (sparse variant_data → hom-ref filled MT → LD prune / matrix / scores).

Variant filtering uses the **Variant Annotation Table (VAT)** ancestry-specific GVS fields when
``--use-vat-callrate`` is set (default): per-ancestry ``gvs_<gen_anc>_an`` / (2 × N) as call rate on
autosomes/PAR, where N is the per-ancestry diploid sample count from
``ld_methods.ANCESTRY_PRED_OTHER_SAMPLE_COUNTS`` (see ``process_vat_for_ld`` in ``ld_methods``).
Cohort genotypes are restricted to samples present in ``META_HT_PATH`` (filtered by ``--gen-anc``).
VAT schema and field definitions:
https://support.researchallofus.org/hc/en-us/articles/4615256690836-Variant-Annotation-Table

Relevant VAT fields (not transcript-specific): ``gvs_<pop>_ac``, ``gvs_<pop>_an``, ``gvs_<pop>_af``,
``gvs_<pop>_sc`` for computed ancestries afr, amr, eas, eur, mid, sas, oth (plus ``gvs_all_*``).
Rows are variant×transcript; we collapse to one row per (locus, alleles).

Design notes (aligned with v8/analysis/gnomad_ld_v4.py):
- Read VDS, optionally filter chromosomes, filter_samples to analysis cohort (meta HT).
- Work from vds.variant_data; row-level ``pop_freq`` for LD still comes from hl.agg.call_stats on
  the filtered cohort MT after VAT row filter.
- With VAT: keep variants passing ancestry call rate from ``gvs_<pop>_an``; without VAT, use in-VDS
  cohort_callrate from call_stats only.
- Optionally require PASS rows when ``filters`` exists (len(filters)==0).
- Checkpoint variant_data intermediate (gnomAD-style), then unfilter_entries + hom-ref fill.
- Optional ``--adj``: gnomAD-style filter **after** hom-ref fill; imputed hom-ref gaps
  (hom-ref ``GT`` with missing ``GQ``) stay passable because they have no stored quality data.

This is **step 1 of 2**. It writes the LD-ready cohort MatrixTable; LD prune / matrix / scores are
step 2, in ``run_ld_from_aou_ld_mt.py``. Both derive their default MT path from ``MY_BUCKET``, so
those two constants must stay in sync. Build the VAT first with ``generate_vat_for_ld.py``.

Example:
  python generate_ld_data_from_aou_vds.py --batch --gen-anc eas --contig chr22 \\
    --generate-ld-mt --min-call-rate 0.9 --vat-path gs://bucket/vat/aou_v8_vat_table.ht \\
    --output-prefix gs://bucket/ld/eas_chr22 --overwrite
"""
from __future__ import annotations

import argparse
import logging
import sys
from typing import Optional

import hail as hl
from hail.utils.misc import new_temp_file

from v8.analysis.ld_methods import (
    ANCESTRY_PRED_OTHER_SAMPLE_COUNTS,
    get_adj_expr as get_adj_expr_aou,
    get_aou_variant_data,
    gen_anc_to_vat_gvs_prefix,
    process_vat_for_ld,
)

logger = logging.getLogger("aou_ld_from_vds")


def configure_logging(level: int = logging.INFO) -> None:
    """Emit INFO logs to stdout so they land in Dataproc/Batch driver output.

    Without this, ``logger.setLevel`` alone leaves the root logger without a handler and every
    progress message in this script is silently dropped.
    """
    logging.basicConfig(
        level=level,
        format="%(asctime)s %(levelname)s %(name)s - %(message)s",
        stream=sys.stdout,
        force=True,
    )
    logger.setLevel(level)

# ---------------------------------------------------------------------------
# Paths and LD frequency defaults (gnomAD LD resources style)
# ---------------------------------------------------------------------------

GEN_ANC_PATH = "gs://fc-aou-datasets-controlled/v8/wgs/short_read/snpindel/aux/ancestry/ancestry_preds.tsv"
VDS_PATH = "gs://fc-aou-datasets-controlled/v8/wgs/short_read/snpindel/vds/hail.vds"
META_HT_PATH = "gs://aou_analysis/v8/data/utils/aou_v8_sample_meta.ht"

VAT_HT_PATH = "gs://aou_marten/ld_data/vat/aou_v8_vat_table.ht"


# Hail scratch (short retention is fine; these are throwaway intermediates).
TMP_BUCKET = "gs://aou_tmp/ld_marten"
# Pipeline outputs. MUST match MY_BUCKET in run_ld_from_aou_ld_mt.py, which derives its default
# --ld-mt-path from it: a mismatch means step 2 looks for an MT step 1 never wrote there.
MY_BUCKET = "gs://aou_tmp_30_days/ld_marten"

def annotate_adj_ld(
    mt: hl.MatrixTable,
    adj_gq: int = 30,
    adj_ab: float = 0.2,
) -> hl.MatrixTable:
    """gnomAD-style ``adj`` for LD after :func:`fill_missing_gt_hom_ref`.

    Hom-ref cells with **missing** ``GQ`` pass ``adj`` in addition to the usual ``get_adj_expr`` rule.
    Those are the sparse matrix gaps filled with hom-ref, which have no stored quality scores; the
    strict rule (GQ ≥ threshold) would drop them all. Real hom-ref calls with non-missing low ``GQ``
    still fail the strict branch.
    """
    if "GT" not in mt.entry and "LGT" in mt.entry:
        gt_expr = mt.LGT
    else:
        gt_expr = mt.GT
    if "AD" not in mt.entry and "LAD" in mt.entry:
        ad_expr = mt.LAD
    else:
        ad_expr = mt.AD
    strict = get_adj_expr_aou(gt_expr, mt.GQ, ad_expr, adj_gq, adj_ab)
    imputed_hom_ref_gap = gt_expr.is_hom_ref() & hl.is_missing(mt.GQ)
    # Missing strict (e.g. het with missing GQ) must be False, not missing, or filter_entries drops it.
    return mt.annotate_entries(adj=hl.coalesce(strict, False) | imputed_hom_ref_gap)


def fill_missing_gt_hom_ref(mt: hl.MatrixTable) -> hl.MatrixTable:
    """Pseudo-dense-ify genotype calls for LD: restore filtered (sparse) entries, then set missing GT to hom-ref.

    VDS variant_data is sparse: many sample×variant pairs have no stored call. ``unfilter_entries``
    brings those cells back into the matrix; any still-missing ``GT`` is filled with ``0/0`` so downstream
    LD (correlations, pruning) sees a complete diploid call matrix.

    ASSUMPTION: restored calls are hom-ref, or enough are hom-ref to be confident in our LD evaluation.

    Args:
        mt: Matrix table whose entries include ``GT`` (typically AoU ``variant_data`` after cohort filters).

    Returns:
        Same schema with every entry row present; ``GT`` is never missing (hom-ref where it was absent).
    """
    logger.info("Starting fill_missing_gt_hom_ref")
    mt = mt.unfilter_entries()
    return mt.annotate_entries(
        GT=hl.if_else(
            hl.is_missing(mt.GT),
            hl.if_else(mt.locus.in_autosome_or_par(), hl.call(0, 0), hl.call(0)),
            mt.GT,
        )
    )


def generate_ld_mt_from_vds(
    vds_path: str,
    meta_ht: hl.Table,
    gen_anc: str,
    contig: Optional[str],
    min_call_rate: float,
    ld_mt_path: str,
    overwrite: bool,
    adj: bool,
    naive_coalesce: Optional[int],
    vat_table: Optional[hl.Table],
    recalculate_callrate_on_other: bool = False,
) -> hl.MatrixTable:
    """
    gnomAD-style: read VDS → filter samples → checkpoint variant_data → row filters
    (VAT ancestry ``gvs_*_an`` call rate and/or in-VDS call_stats) → hom-ref fill
    → optional LD ``adj`` (:func:`annotate_adj_ld`) → checkpoint LD MT.

    ``adj`` runs **after** fill so every cell has ``GT``; :func:`annotate_adj_ld` extends the usual
    gnomAD rule so hom-ref with **missing** ``GQ`` (imputed sparse gaps) can pass, while real calls
    with defined low ``GQ`` or bad het balance still fail.
    """

    gen_anc_key = gen_anc_to_vat_gvs_prefix(gen_anc)
    gen_anc_threshold_str = (
        f"gvs_{gen_anc_key}_callrate_threshold_{str(min_call_rate)}"
    )
    gen_anc_callrate = f"gvs_{gen_anc_key}_callrate"

    gen_anc_af_str = f"gvs_{gen_anc_key}_af"
    gen_anc_an_str = f"gvs_{gen_anc_key}_an"
    gen_anc_ac_str = f"gvs_{gen_anc_key}_ac"

    if contig:
        vat_table = vat_table.filter(vat_table.locus.contig == contig)
    else:
        logger.info("No --contig provided; using all contigs from VAT and VDS.")

    if recalculate_callrate_on_other:

        logger.info(
            "Recalculating %s call rate using ANCESTRY_PRED_OTHER_SAMPLE_COUNTS denominator",
            gen_anc_key,
        )

        sample_count_with_other_correction = ANCESTRY_PRED_OTHER_SAMPLE_COUNTS[
            gen_anc_key
        ]

        callrate_fields = {
            gen_anc_callrate: hl.if_else(
                vat_table.in_autosome_or_par,
                vat_table[f"gvs_{gen_anc_key}_an"]
                / hl.literal(float(sample_count_with_other_correction * 2)),
                vat_table[f"gvs_{gen_anc_key}_an"]
                / hl.literal(float(sample_count_with_other_correction)),
            )
        }
        vat_table = vat_table.annotate(**callrate_fields)

        callrate_threshold_fields = {
            gen_anc_threshold_str: vat_table[gen_anc_callrate] >= min_call_rate
        }
        vat_table = vat_table.annotate(**callrate_threshold_fields)

    # Fail loudly if the VAT was built with a different --min-call-rate: the threshold value is
    # baked into the field name, so a mismatch would otherwise surface as an opaque Hail
    # "no field" error deep in the query.
    if gen_anc_threshold_str not in vat_table.row:
        available = sorted(
            f for f in vat_table.row if "_callrate_threshold_" in f
        )
        raise ValueError(
            f"VAT table has no field {gen_anc_threshold_str!r} (implied by "
            f"--min-call-rate {min_call_rate} and --gen-anc {gen_anc}). "
            f"Threshold fields present: {available or 'none'}. "
            "Rebuild the VAT with generate_vat_for_ld.py --min-call-rate "
            f"{min_call_rate}, or pass --recalculate-callrate-on-other to derive it here."
        )

    vat_table = vat_table.filter(vat_table[f"{gen_anc_threshold_str}"])

    entries_to_keep = ["GT"]

    if adj:
        entries_to_keep.append("GQ")
        entries_to_keep.append("LAD")

    vmt = get_aou_variant_data(
        split=True,
        filter_samples_ht=meta_ht,
        filter_variant_ht=vat_table,
        chrom=[contig] if contig else None,
        naive_coalesce_partitions=naive_coalesce,
        entries_to_keep=entries_to_keep,
        checkpoint_variant_data=False,
        vds_path=vds_path,
    )

    logger.info("checkpoint variant_data (gnomAD-style temp)")

    count_outs = vmt.count()
    if count_outs[0] == 0:
        raise ValueError(
            f"No variants in VDS after VAT filter for gen_anc={gen_anc} and contig={contig}"
        )
    if count_outs[1] == 0:
        raise ValueError(
            f"No samples in VDS after meta filter for gen_anc={gen_anc} and contig={contig}"
        )

    vmt = vmt.checkpoint(new_temp_file("aou_vds_variant", "mt"))

    vmt = vmt.annotate_rows(
        gen_anc_af=vat_table[vmt.row_key][f"{gen_anc_af_str}"],
        gen_anc_an=vat_table[vmt.row_key][f"{gen_anc_an_str}"],
        gen_anc_ac=vat_table[vmt.row_key][f"{gen_anc_ac_str}"],
        gen_anc_callrate=vat_table[vmt.row_key][f"{gen_anc_callrate}"],
    )

    # No need to filter rows based on call rate since we are using the filtered VAT table

    logger.info("unfilter_entries + hom-ref fill → GT")
    vmt = fill_missing_gt_hom_ref(vmt)

    if adj:
        logger.info("annotate_adj_ld + filter_entries (after hom-ref fill)")
        vmt = annotate_adj_ld(vmt)
        vmt = vmt.filter_entries(vmt.adj)
        vmt = vmt.select_entries("GT")

    logger.info("checkpoint LD MT intermediate (pre-callstats)")
    # vmt = vmt.checkpoint(new_temp_file("aou_ld_mt.pre_callstats", "mt"))

    logger.info("Recomputing cohort call_stats into row field gen_anc_freq")
    vmt = vmt.annotate_rows(gen_anc_freq=hl.agg.call_stats(vmt.GT, vmt.alleles))

    logger.info("checkpoint LD MT → %s", ld_mt_path)
    vmt = vmt.checkpoint(ld_mt_path, overwrite=overwrite)

    return vmt


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
            log="/aou_ld_from_vds.log",
            app_name="aou_ld_from_vds",
            gcs_requester_pays_configuration=args.requester_pays,
        )
    elif args.dataproc:
        hl.init(
            tmp_dir=args.tmp_dir,
            default_reference="GRCh38",
            app_name="aou_ld_from_vds",
            gcs_requester_pays_configuration=args.requester_pays,
            backend="spark",
        )
    else:
        raise ValueError(
            "Must run in batch or dataproc mode. Please set --batch or --dataproc."
        )

    hl._set_flags(use_new_shuffle="1")

    # Read (or build) the Variant Annotation Table used to filter to LD-eligible variants.
    # Prefer building it with generate_vat_for_ld.py and passing --vat-path here; --do-generate-vat
    # rebuilds it inline, which repeats expensive work on every LD run.
    vat_path = args.vat_path or VAT_HT_PATH

    vat_table = None
    if args.do_generate_vat:
        logger.info("Building VAT inline and checkpointing to %s", vat_path)
        vat_table = process_vat_for_ld(
            callrate_threshold=args.min_call_rate,
            contig=contig,
        )
        vat_table = vat_table.checkpoint(vat_path, overwrite=args.overwrite)

    if vat_table is None:
        logger.info("Reading VAT from %s", vat_path)
        vat_table = hl.read_table(vat_path)

    # Read in Wenhan's generated sample meta HT. Used to to filter to samples.
    meta_ht = hl.read_table(META_HT_PATH)
    gen_anc_upper = gen_anc.upper()

    # Wenhan's table is capitalized.
    meta_ht_filter = meta_ht.filter(meta_ht.ancestry == gen_anc_upper)
    if meta_ht_filter.count() == 0:
        raise ValueError(f"No samples in meta for gen_anc={gen_anc}")

    # Note that we do not read in the VDS here.

    # Default output prefix and defining some paths.
    output_prefix = (
        args.output_prefix or f"{MY_BUCKET}/{gen_anc}_{contig or 'genome'}_vdsld"
    )

    custom_suffix = args.custom_suffix or ""
    suffix_token = f"_{custom_suffix}" if custom_suffix else ""
    ld_mt_path = args.ld_mt_path or f"{output_prefix}_ld_cohort{suffix_token}.mt"

    if args.generate_ld_mt:
        generate_ld_mt_from_vds(
            vds_path=args.vds_path or VDS_PATH,
            meta_ht=meta_ht_filter,
            gen_anc=gen_anc,
            contig=contig,
            min_call_rate=args.min_call_rate,
            ld_mt_path=ld_mt_path,
            overwrite=args.overwrite,
            adj=args.adj,
            naive_coalesce=args.naive_coalesce,
            vat_table=vat_table,
            recalculate_callrate_on_other=args.recalculate_callrate_on_other,
        )


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="AoU VDS → cohort LD MatrixTable"
    )
    parser.add_argument("--gen-anc", default="eas", help="Cohort pop / meta_ht.pop")
    parser.add_argument(
        "--contig", default=None, help="Restrict to one chrom (e.g. chr22)"
    )
    parser.add_argument(
        "--vds-path",
        default=None,
        help=f"VDS URI (default: ld_methods.VDS_PATH = {VDS_PATH})",
    )
    parser.add_argument(
        "--output-prefix",
        default=None,
        help="Prefix for pruned/index/bm/ldscores paths and default LD MT name",
    )
    parser.add_argument(
        "--ld-mt-path",
        default=None,
        help="Checkpoint path for cohort LD MatrixTable",
    )
    parser.add_argument(
        "--min-call-rate",
        dest="min_call_rate",
        type=float,
        default=0.90,
        help=(
            "Minimum call rate (0–1) used for VAT callrate-threshold filtering."
        ),
    )
    parser.add_argument(
        "--generate-ld-mt",
        action="store_true",
        help="Read VDS, filter, call-stats, hom-ref fill, checkpoint LD MT",
    )

    parser.add_argument(
        "--adj",
        action="store_true",
        help="Keep ADJ genotypes only (needs GQ/DP/AD entries)",
    )
    parser.add_argument(
        "--naive-coalesce",
        type=int,
        default=None,
        help="Coalesce variant_data to N partitions before checkpoint",
    )
    parser.add_argument(
        "--custom-suffix", default="", help="Suffix token in output paths"
    )
    parser.add_argument("--overwrite", action="store_true")
    parser.add_argument("--dataproc", action="store_true")
    parser.add_argument(
        "--batch",
        action="store_true",
        help="Run the pipeline in QueryOnBatch, not Dataproc.",
    )
    parser.add_argument(
        "--test", action="store_true", help="Use chr22 if --contig not set"
    )
    parser.add_argument("--tmp-dir", default=TMP_BUCKET)
    parser.add_argument(
        "--vat-path",
        default=None,
        help=(
            "Processed VAT Hail table to filter variants with, as written by "
            f"generate_vat_for_ld.py (default: {VAT_HT_PATH}). For a single-contig VAT, point "
            "this at e.g. .../aou_v8_vat_table_chr1.ht."
        ),
    )
    parser.add_argument(
        "--do-generate-vat",
        action="store_true",
        help=(
            "Build the VAT inline and checkpoint it to --vat-path before running. Prefer "
            "generate_vat_for_ld.py; this repeats expensive work on every LD run."
        ),
    )
    parser.add_argument(
        "--requester-pays",
        default="aou-neale-gwas",
        help="GCS requester pays project",
    )
    parser.add_argument(
        "--recalculate-callrate-on-other",
        action="store_true",
        help="Recalculate call rate on samples with correction for samples labeled as other in ancestry prediction",
    )

    args = parser.parse_args()

    main(args)
