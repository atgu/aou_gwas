#!/usr/bin/env python3
"""
Generate the AoU Variant Annotation Table (VAT) Hail table used to filter variants for LD.

Supersedes ``generate_vat_chr1.py``, which was hardcoded to one chromosome. The processing
itself is unchanged and lives in :func:`v8.analysis.ld_methods.process_vat_for_ld`: import the
controlled-release VAT TSV, cast the string numeric columns, key by ``(locus, alleles)``,
collapse the variant x transcript duplicates with ``distinct()``, then annotate per-ancestry
``gvs_<anc>_callrate`` and ``gvs_<anc>_callrate_threshold_<value>`` flags.

Two modes:

``--mode single``
    One ``import_table`` pass over the whole VAT, one output table. Fewer, larger jobs.

``--mode serial``
    Loop contigs, writing one table per contig. Each contig is an independent write, so a
    failure costs one contig instead of the whole genome, and completed contigs are skipped on
    re-run (see ``--reuse-existing``). Add ``--union`` to also concatenate the per-contig
    tables into a single genome-wide table at the end.

Which mode to use: the expensive step is ``distinct()``, which shuffles the full
variant x transcript row set. Serial mode shrinks that shuffle to roughly one chromosome's
worth of rows, which is the difference between fitting in worker memory and not. The cost is
that ``hl.import_table`` re-scans and re-decompresses the entire single-file bgz TSV once per
contig -- serial mode does N times the import I/O of single mode. Prefer ``single`` if it
completes; fall back to ``serial`` when the genome-wide shuffle is the thing that fails.

NOTE ON FIELD NAMES: the threshold flag embeds the float, e.g. ``gvs_eur_callrate_threshold_0.9``.
The ``.`` means Hail attribute access does not work -- use ``ht['gvs_eur_callrate_threshold_0.9']``,
not ``ht.gvs_eur_callrate_threshold_0.9``. Keep ``--min-call-rate`` consistent with the value
passed to ``generate_ld_data_from_aou_vds.py``, which reconstructs this field name to filter rows.

VAT schema:
https://support.researchallofus.org/hc/en-us/articles/4615256690836-Variant-Annotation-Table

Examples:
  # One genome-wide table (Query on Batch)
  python generate_vat_for_ld.py --batch --mode single --overwrite

  # Per-contig tables, resumable, plus a unioned genome-wide table
  python generate_vat_for_ld.py --dataproc --mode serial --union --overwrite

  # Just chr1, matching the legacy generate_vat_chr1.py output path
  python generate_vat_for_ld.py --batch --mode serial --contigs chr1 --overwrite
"""
from __future__ import annotations

import argparse
import logging
import sys
from typing import List, Optional

import hail as hl

from v8.analysis.ld_methods import (
    ANCESTRY_PRED_OTHER_SAMPLE_COUNTS,
    VAT_PATH,
    process_vat_for_ld,
)

logger = logging.getLogger("generate_vat_for_ld")

# ---------------------------------------------------------------------------
# Paths and defaults
# ---------------------------------------------------------------------------

# Output root and basename chosen so per-contig output matches the table already written by the
# legacy generate_vat_chr1.py: gs://aou_marten/ld_data/vat/aou_v8_vat_table_chr1.ht
DEFAULT_OUTPUT_ROOT = "gs://aou_marten/ld_data/vat"
DEFAULT_BASENAME = "aou_v8_vat_table"
TMP_BUCKET = "gs://aou_tmp/ld_marten"

# GRCh38 primary contigs. chrM is excluded by default: LD on the mitochondrion is not
# meaningful, and it is not autosome-or-PAR so it would take the haploid callrate denominator.
# Pass it explicitly via --contigs if you want it.
DEFAULT_CONTIGS = [f"chr{i}" for i in range(1, 23)] + ["chrX", "chrY"]


def configure_logging(level: int = logging.INFO) -> None:
    """Emit INFO logs to stdout so they land in Dataproc/Batch driver output."""
    logging.basicConfig(
        level=level,
        format="%(asctime)s %(levelname)s %(name)s - %(message)s",
        stream=sys.stdout,
        force=True,
    )
    logger.setLevel(level)


def vat_ht_path(output_root: str, basename: str, contig: Optional[str] = None) -> str:
    """Output URI for the whole-genome table (``contig=None``) or one contig."""
    suffix = f"_{contig}" if contig else ""
    return f"{output_root}/{basename}{suffix}.ht"


def _write_completed(path: str) -> bool:
    """True if ``path`` looks like a fully written Hail table.

    Checks for the ``_SUCCESS`` marker rather than the directory itself, so a partially written
    table left behind by a killed job is not mistaken for a finished one.
    """
    try:
        return hl.hadoop_exists(f"{path}/_SUCCESS")
    except Exception as e:  # noqa: BLE001 - unreadable path is treated as "not written"
        logger.warning("Could not stat %s (%s); treating as not written", path, e)
        return False


def generate_vat(
    contig: Optional[str],
    output_path: str,
    vat_path: str,
    callrate_threshold: float,
    overwrite: bool,
) -> hl.Table:
    """Process the VAT (optionally for one contig) and write it to ``output_path``."""
    logger.info(
        "Processing VAT for %s -> %s", contig or "all contigs", output_path
    )
    vat_table = process_vat_for_ld(
        n_samples_dict=ANCESTRY_PRED_OTHER_SAMPLE_COUNTS,
        vat_path=vat_path,
        callrate_threshold=callrate_threshold,
        contig=contig,
    )
    vat_table = vat_table.checkpoint(output_path, overwrite=overwrite)
    logger.info("Wrote %s (%s rows)", output_path, vat_table.count())
    return vat_table


def run_single(args: argparse.Namespace) -> None:
    """One pass over the whole VAT, one genome-wide output table."""
    out = args.output or vat_ht_path(args.output_root, args.basename)
    if args.reuse_existing and _write_completed(out):
        logger.info("%s already written; nothing to do (--reuse-existing)", out)
        return
    generate_vat(
        contig=None,
        output_path=out,
        vat_path=args.vat_path,
        callrate_threshold=args.min_call_rate,
        overwrite=args.overwrite,
    )


def run_serial(args: argparse.Namespace, contigs: List[str]) -> None:
    """One output table per contig, then optionally union them into one table."""
    written: List[str] = []
    skipped: List[str] = []
    failed: List[str] = []

    for contig in contigs:
        out = vat_ht_path(args.output_root, args.basename, contig)
        if args.reuse_existing and _write_completed(out):
            logger.info("Skipping %s: %s already written", contig, out)
            skipped.append(contig)
            written.append(out)
            continue
        try:
            generate_vat(
                contig=contig,
                output_path=out,
                vat_path=args.vat_path,
                callrate_threshold=args.min_call_rate,
                overwrite=args.overwrite,
            )
            written.append(out)
        except Exception as e:  # noqa: BLE001 - keep going so one bad contig is not fatal
            if not args.continue_on_error:
                raise
            logger.error("Contig %s FAILED: %s", contig, e)
            failed.append(contig)

    logger.info(
        "Serial run done: %s written, %s reused, %s failed",
        len(written) - len(skipped),
        len(skipped),
        len(failed),
    )
    if failed:
        logger.error("Failed contigs (rerun with --reuse-existing to fill gaps): %s", failed)

    if not args.union:
        return
    if failed:
        raise ValueError(
            f"Refusing to --union an incomplete set; these contigs failed: {failed}. "
            "Rerun with --reuse-existing to fill the gaps, then --union."
        )

    union_out = args.output or vat_ht_path(args.output_root, args.basename)
    logger.info("Unioning %s per-contig tables into %s", len(written), union_out)
    # Contigs are processed in reference order and each table is keyed by (locus, alleles),
    # so this concatenation is already key-sorted.
    tables = [hl.read_table(p) for p in written]
    unioned = tables[0] if len(tables) == 1 else tables[0].union(*tables[1:])
    unioned = unioned.checkpoint(union_out, overwrite=args.overwrite)
    logger.info("Wrote %s (%s rows)", union_out, unioned.count())


def main(args: argparse.Namespace) -> None:
    configure_logging()

    contigs = args.contigs or DEFAULT_CONTIGS
    if args.mode == "single" and args.contigs:
        raise ValueError(
            "--contigs only applies to --mode serial. For a contig subset in one table, run "
            "--mode serial --contigs ... --union."
        )
    if args.union and args.mode != "serial":
        raise ValueError("--union only applies to --mode serial.")

    if args.batch:
        hl.init(
            backend="batch",
            tmp_dir=args.tmp_dir,
            driver_memory="highmem",
            driver_cores=8,
            worker_memory="standard",
            worker_cores=8,
            default_reference="GRCh38",
            log="/generate_vat_for_ld.log",
            app_name="generate_vat_for_ld",
            gcs_requester_pays_configuration=args.requester_pays,
        )
    elif args.dataproc:
        hl.init(
            tmp_dir=args.tmp_dir,
            default_reference="GRCh38",
            app_name="generate_vat_for_ld",
            gcs_requester_pays_configuration=args.requester_pays,
            backend="spark",
        )
    else:
        raise ValueError("Must run in batch or dataproc mode. Set --batch or --dataproc.")

    hl._set_flags(use_new_shuffle="1")

    if args.mode == "single":
        run_single(args)
    else:
        logger.info("Serial mode over %s contigs: %s", len(contigs), contigs)
        run_serial(args, contigs)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description=(
            "Generate the AoU VAT Hail table for LD variant filtering, either as one "
            "genome-wide table (--mode single) or one table per contig (--mode serial)."
        )
    )
    parser.add_argument(
        "--mode",
        choices=["single", "serial"],
        default="single",
        help=(
            "single: one pass over the whole VAT into one table. "
            "serial: one table per contig, resumable, optionally unioned with --union. "
            "(default: single)"
        ),
    )
    parser.add_argument(
        "--contigs",
        nargs="+",
        default=None,
        metavar="CONTIG",
        help=(
            "Contigs for --mode serial (e.g. --contigs chr1 chr2 chrX). "
            f"Default: {' '.join(DEFAULT_CONTIGS)}"
        ),
    )
    parser.add_argument(
        "--union",
        action="store_true",
        help=(
            "After --mode serial, concatenate the per-contig tables into one genome-wide "
            "table at --output (or the default whole-genome path)."
        ),
    )
    parser.add_argument(
        "--vat-path",
        default=VAT_PATH,
        help=f"URI of the AoU VAT TSV (default: {VAT_PATH})",
    )
    parser.add_argument(
        "--output-root",
        default=DEFAULT_OUTPUT_ROOT,
        help=f"Directory for outputs (default: {DEFAULT_OUTPUT_ROOT})",
    )
    parser.add_argument(
        "--basename",
        default=DEFAULT_BASENAME,
        help=(
            "Output table basename; per-contig tables get a _<contig> suffix "
            f"(default: {DEFAULT_BASENAME})"
        ),
    )
    parser.add_argument(
        "--output",
        default=None,
        help=(
            "Explicit path for the single/unioned genome-wide table. Overrides "
            "--output-root/--basename. Does not affect per-contig paths."
        ),
    )
    parser.add_argument(
        "--min-call-rate",
        dest="min_call_rate",
        type=float,
        default=0.90,
        help=(
            "Call-rate threshold baked into the gvs_<anc>_callrate_threshold_<value> field "
            "names. Must match the --min-call-rate passed to generate_ld_data_from_aou_vds.py. "
            "(default: 0.90 -> field suffix '0.9')"
        ),
    )
    parser.add_argument(
        "--reuse-existing",
        action="store_true",
        help=(
            "Skip any output whose _SUCCESS marker already exists. Use this to resume an "
            "interrupted serial run without recomputing finished contigs."
        ),
    )
    parser.add_argument(
        "--continue-on-error",
        action="store_true",
        help=(
            "In serial mode, log and continue past a failing contig instead of aborting. "
            "--union will still refuse to run if any contig failed."
        ),
    )
    parser.add_argument("--overwrite", action="store_true")
    parser.add_argument("--dataproc", action="store_true")
    parser.add_argument(
        "--batch",
        action="store_true",
        help="Run in Query on Batch (recommended)",
    )
    parser.add_argument("--tmp-dir", default=TMP_BUCKET)
    parser.add_argument(
        "--requester-pays",
        default="aou-neale-gwas",
        help="GCS requester-pays project",
    )

    main(parser.parse_args())
