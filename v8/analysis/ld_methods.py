"""
All of Us LD helpers: VAT (GVS) tables for ancestry filters, and VDS ``variant_data`` loading
with optional multiallelic splitting via ``hl.experimental.sparse_split_multi``.
"""

from __future__ import annotations

import logging
from typing import Optional, List, Union, Set

import hail as hl
from hail.utils.misc import new_temp_file


logger = logging.getLogger("aou_ld_methods")
logger.setLevel(logging.INFO)

VDS_PATH = "gs://fc-aou-datasets-controlled/v8/wgs/short_read/snpindel/vds/hail.vds"


def get_adj_expr(
    gt_expr: hl.expr.CallExpression,
    gq_expr: Union[hl.expr.Int32Expression, hl.expr.Int64Expression],
    ad_expr: hl.expr.ArrayNumericExpression,
    adj_gq: int = 30,
    adj_ab: float = 0.2,
) -> hl.expr.BooleanExpression:
    """
    Get adj genotype annotation. Defaults correspond to gnomAD values.

    Verbatim copy of ``get_adj_expr`` in ``v8/data_prep/process_data.py`` (author: wlu).
    Duplicated deliberately: importing that module executes ``from gnomad.utils.vep import *``,
    ``hailtop.batch``, ``pandas`` and ``tqdm`` at load time, so the LD pipeline would need the
    full GWAS dependency set just to build this one expression. Keep in sync if the upstream
    definition changes.
    """
    return (gq_expr >= adj_gq) & (
        hl.case()
        .when(~gt_expr.is_het(), True)
        .when(gt_expr.is_het_ref(), ad_expr[gt_expr[1]] / hl.sum(ad_expr) >= adj_ab)
        .default(
            (ad_expr[gt_expr[0]] / hl.sum(ad_expr) >= adj_ab)
            & (ad_expr[gt_expr[1]] / hl.sum(ad_expr) >= adj_ab)
        )
    )

AUX_PATH = "gs://fc-aou-datasets-controlled/v8/wgs/short_read/snpindel/aux"
VAT_PATH = f'{AUX_PATH}/vat/vat_complete.bgz.tsv.gz'

# We are not calculating LD on OTH. 
VAT_GVS_ANCESTRY_CODES = frozenset({"afr", "amr", "eas", "eur", "mid", "sas"})

META_GEN_ANC_SAMPLE_COUNTS = {'afr': 77444,
    'amr': 71540,
    'eas': 9488,
    'eur': 227273,
    'mid': 1153,
    'sas': 5132}

ANCESTRY_PRED_OTHER_SAMPLE_COUNTS = {'afr': 79826,
 'amr': 71854,
 'eas': 9440,
 'eur': 223350,
 'mid': 810,
 'oth': 25504,
 'sas': 4046}


def gen_anc_to_vat_gvs_prefix(gen_anc: str) -> str:
    """
    Map a genetic-ancestry code to the VAT ``gvs_<pop>_*`` column prefix.

    :param gen_anc: Short code (e.g. ``eas``, ``afr``) matching ``VAT_GVS_ANCESTRY_CODES``.
    :return: Lowercased ancestry token used in column names (e.g. ``eas`` → ``gvs_eas_an``).
    :raises ValueError: If ``gen_anc`` is not a supported GVS ancestry.
    """
    g = gen_anc.strip().lower()
    if g not in VAT_GVS_ANCESTRY_CODES:
        raise ValueError(
            f"gen_anc={gen_anc!r} has no GVS VAT fields; expected one of {sorted(VAT_GVS_ANCESTRY_CODES)}"
        )
    return g


def _vat_str_to_int32(s: hl.expr.StringExpression) -> hl.expr.Int32Expression:
    """Parse VAT TSV integer fields when ``import_table`` leaves them as strings."""
    blank = hl.is_missing(s) | (s == "") | (s == ".") | (s == "NA")
    return hl.if_else(blank, hl.missing(hl.tint32), hl.int32(s))


def _vat_str_to_float64(s: hl.expr.StringExpression) -> hl.expr.Float64Expression:
    """Parse VAT TSV float fields (e.g. AF) when stored as strings."""
    blank = hl.is_missing(s) | (s == "") | (s == ".") | (s == "NA")
    return hl.if_else(blank, hl.missing(hl.tfloat64), hl.float64(s))


def process_vat_for_ld(
    n_samples_dict: dict[str, int] = ANCESTRY_PRED_OTHER_SAMPLE_COUNTS,
    vat_path: str = VAT_PATH,
    callrate_threshold: float = 0.90,
    test: bool = False,
    contig: Optional[str] = None,
) -> hl.Table:
    """
    Process an AoU Variant Annotation Table (Hail) for LD analysis.

    :param vat_path: URI of the VAT TSV (``import_table`` without ``impute=True`` leaves
        numeric columns as strings; those are cast to int/float before call-rate math).
    :param n_samples_dict: Per-ancestry diploid sample count ``N`` per ancestry for
        ``gvs_<anc>_callrate = gvs_<anc>_an / (2 * N)`` on autosomes/PAR (haploid sex
        chromosomes use ``/ N``). Defaults to :data:`ANCESTRY_PRED_OTHER_SAMPLE_COUNTS`
        (denominators aligned with ancestry prediction / “other” cohort sizes), not
        :data:`META_GEN_ANC_SAMPLE_COUNTS`.
    :param callrate_threshold: Minimum call rate threshold for variants. Becomes part of the
        output field name: ``gvs_<anc>_callrate_threshold_<callrate_threshold>``.
    :param test: Shorthand for ``contig="chr1"``. Ignored when ``contig`` is set.
    :param contig: Keep only this contig (e.g. ``chr1``). ``None`` processes the whole VAT.
        The filter is applied on the raw ``contig`` string column, before ``hl.locus`` is
        built, so it costs nothing beyond the TSV scan.
    :return: Keyed table joinable to ``variant_data`` on ``(locus, alleles)``.
    """
    vat_table = hl.import_table(vat_path, quote='"', delimiter="\t", force_bgz=True)

    if contig is None and test:
        contig = "chr1"
    if contig is not None:
        logger.info("Filtering VAT to contig %s", contig)
        vat_table = vat_table.filter(vat_table.contig == contig)

    vat_table = vat_table.annotate(
        locus=hl.locus(vat_table.contig, _vat_str_to_int32(vat_table.position)),
        alleles=hl.array([vat_table.ref_allele, vat_table.alt_allele]),
    )

    gvs_metric = ["ac", "an", "af", "sc"]
    gvs_genanc = list(VAT_GVS_ANCESTRY_CODES) + ["max", "all"]

    gvs_rows = [
        f"gvs_{gvs_genanc_i}_{gvs_metric_i}"
        for gvs_genanc_i in gvs_genanc
        for gvs_metric_i in gvs_metric
    ]

    filter_rows = [
        "locus",
        "alleles",
        "contig",
        "position",
        "ref_allele",
        "alt_allele",
    ] + gvs_rows

    vat_table = vat_table.select(*filter_rows)

    # ``import_table`` without ``impute=True`` leaves numeric VAT columns as strings; cast before math.
    vat_table = vat_table.select(
        vat_table.locus,
        vat_table.alleles,
        vat_table.contig,
        vat_table.position,
        vat_table.ref_allele,
        vat_table.alt_allele,
        **{
            f"gvs_{g}_{m}": (
                _vat_str_to_float64(vat_table[f"gvs_{g}_{m}"])
                if m == "af"
                else _vat_str_to_int32(vat_table[f"gvs_{g}_{m}"])
            )
            for g in gvs_genanc
            for m in gvs_metric
        },
    )
    
    vat_table = vat_table.checkpoint(new_temp_file('locus_allele_vat',extension='ht'))
    
    vat_table = vat_table.key_by('locus','alleles')

    # There is no difference between callstats for the same variant across different transcripts, so we can call distinct

    vat_table = vat_table.distinct()

    vat_table = vat_table.annotate(
        in_autosome_or_par = vat_table.locus.in_autosome_or_par()
    )

    callrate_fields = {
        f"gvs_{gen_anc}_callrate": hl.if_else(
            vat_table.in_autosome_or_par,
            vat_table[f"gvs_{gen_anc}_an"] / hl.literal(float(n_samples_dict[gen_anc] * 2)),
            vat_table[f"gvs_{gen_anc}_an"] / hl.literal(float(n_samples_dict[gen_anc])),  # haploid sex chromosomes
        )
        for gen_anc in VAT_GVS_ANCESTRY_CODES
    }

    vat_table = vat_table.annotate(**callrate_fields)

    callrate_threshold_fields = {
        f"gvs_{gen_anc}_callrate_threshold_{str(callrate_threshold)}": vat_table[f"gvs_{gen_anc}_callrate"] >= callrate_threshold
        for gen_anc in VAT_GVS_ANCESTRY_CODES
    }
    vat_table = vat_table.annotate(**callrate_threshold_fields)

    return vat_table


def _split_and_filter_variant_data_for_loading(
    mt: hl.MatrixTable,
    filter_variant_ht: Optional[hl.Table] = None,
    entries_to_keep: Optional[List[str]] = None,
    checkpoint_before_split: bool = False,
) -> hl.MatrixTable:
    """
    Prepare sparse ``variant_data`` for splitting: optional entry subset, locus pre-filter,
    row annotations, checkpoint, then ``hl.experimental.sparse_split_multi``.

    Order of operations:

    1. If ``entries_to_keep`` is set, select those fields (mapping ``GT``/``AD``/``PL`` to
       ``LGT``/``LAD``/``LPL`` and always keeping ``LA``) before split.
    2. If ``filter_variant_ht`` is set, restrict rows by **locus** only (sorted key-by
       optimization), then annotate ``n_unsplit_alleles`` and ``mixed_site``.
    3. Optional checkpoint (materialize before split).
    4. ``sparse_split_multi(..., filter_changed_loci=True)``.
    5. If ``filter_variant_ht`` is set, ``semi_join_rows`` to the full variant key table.

    :param mt: Unsplit ``vds.variant_data`` MatrixTable.
    :param filter_variant_ht: Optional table keyed by at least ``locus`` (first pass) and
        by ``(locus, alleles)`` after split for ``semi_join_rows``.
    :param entries_to_keep: Names of **global** entry fields to retain before split
        (e.g. ``GT``); these are converted to local ``L*`` where required for splitting.
    :param checkpoint_before_split: If True, checkpoint immediately before ``sparse_split_multi``.
    :return: Split, biallelic-oriented ``variant_data`` MatrixTable.
    """
    if entries_to_keep is not None:
        split_entries_to_keep = entries_to_keep + ["LA"]
        split_entries_to_keep = [
            "L" + e if e in {"GT", "AD", "PL"} else e for e in split_entries_to_keep
        ]
        mt = mt.select_entries(*split_entries_to_keep)

    if filter_variant_ht is not None:
        # Prevents hail from running sort on HT which is already sorted.
        filter_locus_ht = hl.Table(
            hl.ir.TableKeyBy(filter_variant_ht._tir, ["locus"], is_sorted=True)
        )
        mt = mt.filter_rows(hl.is_defined(filter_locus_ht[mt.locus]))

    mt = mt.annotate_rows(
        n_unsplit_alleles=hl.len(mt.alleles),
        mixed_site=(hl.len(mt.alleles) > 2)
        & hl.any(lambda a: hl.is_indel(mt.alleles[0], a), mt.alleles[1:])
        & hl.any(lambda a: hl.is_snp(mt.alleles[0], a), mt.alleles[1:]),
    )

    if checkpoint_before_split:
        mt = mt.checkpoint(new_temp_file("variant_data.before_split", "mt"))

    logger.info("Splitting multiallelics...")
    mt = hl.experimental.sparse_split_multi(mt, filter_changed_loci=True)

    # Different considerations before splitting and after splitting.
    # Note that this does not annotate relevant fields needed for later. 
    if filter_variant_ht is not None:
        mt = mt.semi_join_rows(filter_variant_ht)

    return mt

def get_aou_variant_data(
    split: bool = False,
    filter_samples_ht: Optional[hl.Table] = None,
    n_partitions: Optional[int] = None,
    filter_partitions: Optional[List[int]] = None,
    chrom: Optional[Union[str, List[str], Set[str]]] = None,
    filter_variant_ht: Optional[hl.Table] = None,
    filter_intervals: Optional[List[Union[str, hl.tinterval]]] = None,
    split_reference_blocks: bool = True,
    entries_to_keep: Optional[List[str]] = None,
    checkpoint_variant_data: bool = False,
    naive_coalesce_partitions: Optional[int] = None,
    vds_path: str = VDS_PATH,
) -> hl.MatrixTable:
    """
    Load the AoU short-read WGS VDS and return the **variant_data** MatrixTable (sparse).

    This function does **not** densify and does not return reference blocks; downstream
    code should treat the result as sparse ``variant_data`` only.

    :param split: If True, run :func:`_split_and_filter_variant_data_for_loading` (multiallelic split).
    :param filter_samples_ht: Optional sample table passed to ``hl.vds.filter_samples`` (must match VDS column keys).
    :param n_partitions: Subset partitions when reading (see ``n_partitions`` + ``chrom`` branch).
    :param filter_partitions: Restrict to listed partition indices on both reference and variant MTs.
    :param chrom: Contig(s) to keep (string or list); combined with ``autosomes_only`` / ``sex_chr_only``.
    :param autosomes_only: Set ``chrom`` to all autosomes (GRCh38 reference genome).
    :param sex_chr_only: Set ``chrom`` to X/Y contigs only.
    :param filter_variant_ht: When ``split`` is True, forwarded to :func:`_split_and_filter_variant_data_for_loading`.
        Requires ``split=True`` (raises otherwise).
    :param filter_intervals: Genomic intervals for ``hl.vds.filter_intervals``.
    :param split_reference_blocks: Passed through to ``filter_intervals``.
    :param entries_to_keep: Global entry field names; if ``split`` is True, applied inside the split helper;
        if ``split`` is False, ``select_entries`` is applied after loading.
    :param checkpoint_variant_data: If True, passed as ``checkpoint_before_split`` to the split helper when
        ``split`` is True; if True after that, checkpoints again after optional ``select_entries``.
    :param naive_coalesce_partitions: ``naive_coalesce`` both reference and variant MTs before extraction.
    :param vds_path: URI of the VDS to read (defaults to :data:`VDS_PATH`).
    :return: The **variant_data** ``MatrixTable`` (optionally split and filtered).
    :raises ValueError: If ``filter_variant_ht`` is set and ``split`` is False.
    """

    if filter_variant_ht is not None and split is False:
        raise ValueError(
            "Filtering to a specific set of variants is only supported when splitting"
            " the VDS."
        )

    if isinstance(chrom, str):
        chrom = [chrom]

    # This is not a use case for us actually, since we are only doing single chromosomes at a time. 
    # Removing autosomes_only and sex_chr_only. 
    # This is not a use case for us actually, since we are only doing single chromosomes at a time. 
    # Removing autosomes_only and sex_chr_only. 

    if n_partitions and chrom:
        logger.info(
            "Filtering to chromosome(s) %s with %s partitions...", chrom, n_partitions
        )
        reference_data = hl.read_matrix_table(
            hl.vds.VariantDataset._reference_path(vds_path)
        )
        reference_data = hl.filter_intervals(
            reference_data,
            [hl.parse_locus_interval(x, reference_genome="GRCh38") for x in chrom],
        )
        intervals = reference_data._calculate_new_partitions(n_partitions)
        reference_data = hl.read_matrix_table(
            hl.vds.VariantDataset._reference_path(vds_path),
            _intervals=intervals,
        )
        variant_data = hl.read_matrix_table(
            hl.vds.VariantDataset._variants_path(vds_path),
            _intervals=intervals,
        )
        vds = hl.vds.VariantDataset(reference_data, variant_data)
    elif n_partitions:
        vds = hl.vds.read_vds(vds_path, n_partitions=n_partitions)
    else:
        vds = hl.vds.read_vds(vds_path)
        if chrom:
            logger.info("Filtering to chromosome %s...", chrom)
            vds = hl.vds.filter_chromosomes(vds,keep=chrom)

    # Do note that this comes after filtering to chromosomes. 
    if naive_coalesce_partitions:
        vds = hl.vds.VariantDataset(
            vds.reference_data.naive_coalesce(naive_coalesce_partitions),
            vds.variant_data.naive_coalesce(naive_coalesce_partitions),
        )

    if filter_partitions:
        logger.info("Filtering to %s partitions...", len(filter_partitions))
        vds = hl.vds.VariantDataset(
            vds.reference_data._filter_partitions(filter_partitions),
            vds.variant_data._filter_partitions(filter_partitions),
        )

    if filter_intervals:
        logger.info("Filtering to %s intervals...", len(filter_intervals))
        if isinstance(filter_intervals[0], str):
            filter_intervals = [
                hl.parse_locus_interval(x, reference_genome="GRCh38")
                for x in filter_intervals
            ]
        vds = hl.vds.filter_intervals(
            vds, filter_intervals, split_reference_blocks=split_reference_blocks
        )

    # Do NOT explicitly worry about chr19:5787204. Unclear if it is a problem in AOU VDS.

    # Filter to specific samples. Stress on this is done upstream and outside of this function.
    if filter_samples_ht:
        vds = hl.vds.filter_samples(vds, samples=filter_samples_ht,keep=True)


    # Horrible plan: 
    # n_samples = vmt.count_cols()
    # filter_variant_ht = filter_variant_ht.annotate_rows(
        #  recalculated_callrate = filter_variant_ht[f"filter_variant_ht.gvs_{gen_anc.lower()}_an"] / n_samples
    # )

    # From this point on, we are working with the variant data MT. 
    # Since we will NOT densify, we do not need to worry about splitting or filtering the reference data.
    vmt = vds.variant_data

    # Checkpoint before this split is DIFFERENT than checkpointing after this is split. :brain 
    if split:
        vmt = _split_and_filter_variant_data_for_loading(
            vmt, filter_variant_ht, entries_to_keep, checkpoint_before_split=True
        )

    # Our trusted soruce _split_and_filter_variant_data_for_loading() handles this logic. 
    # if entries_to_keep is not None:
    #     vmt = vmt.select_entries(*entries_to_keep)

    if checkpoint_variant_data:
        vmt = vmt.checkpoint(new_temp_file("vds_loading.variant_data", "mt"))

    return vmt