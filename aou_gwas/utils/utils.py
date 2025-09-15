import hail as hl
import pandas as pd
from aou_gwas.utils.resources import *
from typing import Union
import hailtop.fs as hfs
import pickle


root = "gs://aou_wlu"
data_path = f"{root}/data"
test_path = f"{root}/test"



def parse_empty_missing(v, f):
    new_v = f(v)
    return hl.if_else((hl.len(v) == 0), hl.missing(new_v.dtype), new_v)

def write_pickle_dict(output: str, dict: dict):
    with hfs.open(output, "wb") as f:
        pickle.dump(dict, f)
    f.close()

def read_pickle_dict(dict_path: str):
    dict = pd.read_pickle(dict_path)
    return dict

def irnt(he: hl.expr.Expression, output_loc: str = 'irnt'):
    ht = he._indices.source
    n_rows = ht.aggregate(hl.agg.count_where(hl.is_defined(he)))
    print(n_rows)  # count non-missings
    ht = ht.order_by(he).add_index()
    ht = ht.annotate(**{output_loc: hl.qnorm((ht.idx + 0.5) / n_rows)})
    ht = ht.annotate(**{output_loc: hl.or_missing(~hl.is_nan(ht[output_loc]), ht[output_loc])})
    return ht
    # def irnt(he: hl.expr.Expression, output_loc: str = 'irnt'):
    #     ht = he._indices.source
    #     n_rows = ht.aggregate(hl.agg.count_where(hl.is_defined(ht.x))) # count non-missings
    #     ht = ht.order_by(he).add_index()
    #     return ht.annotate(**{output_loc: hl.qnorm((ht.idx + 0.5) / n_rows)})
    #     # 1) add hl.or_missing(hl.qnorm)
    #     # 2) use approx_cdf for later runs


def get_vat_field_types():
    # For more info: https://support.researchallofus.org/hc/en-us/articles/4615256690836
    tags = (
        "afr",
        "amr",
        "asj",
        "eas",
        "eur",
        "fin",
        "mid",
        "nfr",
        "sas",
        "oth",
        "max",
        "all",
    )
    fields = ("ac", "an", "sc")
    types = {}
    for tag in tags:
        types[f"gvs_{tag}_af"] = hl.tfloat64
        types[f"gnomad_{tag}_af"] = hl.tfloat64
        for field in fields:
            types[f"gvs_{tag}_{field}"] = hl.tint32
            types[f"gnomad_{tag}_{field}"] = hl.tint32
    types["is_canonical_transcript"] = hl.tbool
    types["omim_phenotypes_id"] = hl.tarray(hl.tint32)
    for x in [
        "position",
        "gene_omim_id",
        "splice_ai_acceptor_gain_distance",
        "splice_ai_acceptor_loss_distance",
        "splice_ai_donor_gain_distance",
        "splice_ai_donor_loss_distance",
    ]:
        types[x] = hl.tint32
    for x in [
        "revel",
        "splice_ai_acceptor_gain_score",
        "splice_ai_acceptor_loss_score",
        "splice_ai_donor_gain_score",
        "splice_ai_donor_loss_score",
    ]:
        types[x] = hl.tfloat64
    for x in [
        "omim_phenotypes_name",
        "clinvar_classification",
        "clinvar_phenotype",
        "consequence",
        "dbsnp_rsid",
    ]:
        types[x] = hl.tarray(hl.tstr)
    return types


def group_gene_lens(ht, n, overwrite=False):
    if overwrite:
        gtf = hl.experimental.import_gtf(
            gtf_path, reference_genome="GRCh38", skip_invalid_contigs=True
        )
        gtf = gtf.filter(gtf.feature == "gene")
        gtf = gtf.annotate(
            chrom=gtf.interval.start.contig,
            start_pos=gtf.interval.start.position,
            end_pos=gtf.interval.end.position,
        )
        gtf = gtf.annotate(length=gtf.end_pos - gtf.start_pos + 1)
        gtf = gtf.add_index()
        ht = ht.annotate(group_id=ht.idx // n)
        sub = ht.select(
            "gene_name",
            "gene_id",
            "transcript_id",
            "chrom",
            "start_pos",
            "end_pos",
            "length",
            "idx",
            "group_id",
        )
        sub.write(f"{data_path}/gene_group_by_position_{n}.ht", overwrite=overwrite)
    ht = hl.read_table(f"{data_path}/gene_group_by_position_{n}.ht")
    ht = ht.group_by("group_id").aggregate(group_length=hl.agg.sum(ht.length))
    return ht


def group_gene_interval(gtf, n, path, overwrite=False):
    from collections import Counter
    gtf = gtf.filter((gtf.feature == 'gene') & (gtf.gene_type == 'protein_coding'))
    gtf = gtf.annotate(
        chrom=gtf.interval.start.contig,
        start_pos=gtf.interval.start.position,
        end_pos=gtf.interval.end.position,
    )
    gtf = gtf.annotate(length=gtf.end_pos - gtf.start_pos + 1)
    gtf = gtf.add_index()

    n_gene = gtf.group_by("chrom").aggregate(n_gene_per_chrom=hl.agg.count())
    n_gene = n_gene.annotate(n_gene_per_chrom_round=
                             hl.if_else(n_gene.n_gene_per_chrom % n > 0,
                                        n_gene.n_gene_per_chrom + n - n_gene.n_gene_per_chrom % n,
                                        n_gene.n_gene_per_chrom
                                        ))
    n_gene = n_gene.annotate(chrom_index=n_gene.chrom[3:])
    n_gene = n_gene.annotate(
        chrom_index=hl.case()
        .when(n_gene.chrom_index == "X", "23")
        .when(n_gene.chrom_index == "Y", "24")
        .when(n_gene.chrom_index == "M", "25")
        .default(n_gene.chrom_index)
    )
    n_gene = n_gene.annotate(chrom_index=hl.int64(n_gene.chrom_index))
    n_gene = n_gene.order_by("chrom_index")

    n_gene_lst = n_gene.aggregate(
        hl.cumulative_sum(hl.agg.collect(n_gene.n_gene_per_chrom))
    )
    n_gene_lst.insert(0, 0)
    n_gene_lst.pop()
    print('------- Cumulative number of genes per chromosome: ---------')
    print(n_gene_lst)

    n_gene_round_lst = n_gene.aggregate(hl.cumulative_sum(hl.agg.collect(n_gene.n_gene_per_chrom_round)))
    n_gene_round_lst.insert(0, 0)
    n_gene_round_lst.pop()
    print('------- Cumulative number of genes per chromosome (rounded): ---------')
    print(n_gene_round_lst)


    chrom_lst = [f"chr{i + 1}" for i in range(22)] + ["chrX", "chrY", "chrM"]
    n_gene_df = pd.DataFrame(data={'chrom': chrom_lst, 'cum_n_gene': n_gene_lst, 'cum_n_gene_round':n_gene_round_lst})
    n_gene_ht = hl.Table.from_pandas(n_gene_df, key="chrom")
    print(n_gene_ht.show(26))

    gtf = gtf.annotate(n_previous_genes = n_gene_ht[gtf.chrom].cum_n_gene,
                       n_previous_genes_round = n_gene_ht[gtf.chrom].cum_n_gene_round,)
    gtf = gtf.annotate(
        new_idx=hl.if_else(
            gtf.n_previous_genes > 0, gtf.idx + gtf.n_previous_genes_round - gtf.n_previous_genes, gtf.idx
        )
    )

    gtf = gtf.annotate(group_id=gtf.new_idx // n)
    group_cnt = gtf.group_by('group_id').aggregate(cnt=hl.agg.count(),
                                                   group_length=hl.agg.sum(gtf.length))
    group_dict = gtf.aggregate(hl.agg.counter(gtf.group_id))
    print(Counter(group_dict.values()))

    gtf = gtf.annotate(n_genes_per_group=group_cnt[gtf.group_id].cnt)
    gtf = gtf.annotate(group_id=hl.if_else(gtf.n_genes_per_group < (n/2), gtf.group_id - 1, gtf.group_id))
    group_dict = gtf.aggregate(hl.agg.counter(gtf.group_id))
    print(Counter(group_dict.values()))

    gtf = gtf.select(
        "gene_name",
        "gene_id",
        "transcript_id",
        "chrom",
        "start_pos",
        "end_pos",
        "length",
        "idx",
        "new_idx",
        "group_id",
        "n_previous_genes",
    )

    sub_gtf = gtf.checkpoint(
        f"{path[:-3]}_tmp.ht",
        _read_if_exists=not overwrite,
        overwrite=overwrite,
    )

    # sub_gtf.show()

    group_ht = sub_gtf.group_by("group_id", "chrom").aggregate(
        group_start=hl.agg.min(sub_gtf.start_pos),
        group_end=hl.agg.max(sub_gtf.end_pos),
        genes=hl.agg.collect(sub_gtf.gene_name),
    )
    group_ht = group_ht.annotate(
        interval=hl.locus_interval(
            group_ht.chrom,
            group_ht.group_start,
            group_ht.group_end,
            reference_genome="GRCh38",
            includes_end=True,
            invalid_missing=True,
        )
    )

    group_ht = group_ht.checkpoint(
        path,
        _read_if_exists=not overwrite,
        overwrite=overwrite,
    )

    return group_ht


def gt_to_gp(mt, location: str = "GP"):
    return mt.annotate_entries(
        **{
            location: hl.or_missing(
                hl.is_defined(mt.GT),
                hl.map(
                    lambda i: hl.cond(
                        mt.GT.unphased_diploid_gt_index() == i, 1.0, 0.0
                    ),
                    hl.range(0, hl.triangle(hl.len(mt.alleles))),
                ),
            )
        }
    )


def impute_missing_gp(mt, location: str = "GP", mean_impute: bool = True):
    mt = mt.annotate_entries(_gp=mt[location])
    if mean_impute:
        mt = mt.annotate_rows(
            _mean_gp=hl.agg.array_agg(lambda x: hl.agg.mean(x), mt._gp)
        )
        gp_expr = mt._mean_gp
    else:
        gp_expr = [1.0, 0.0, 0.0]
    return mt.annotate_entries(**{location: hl.or_else(mt._gp, gp_expr)}).drop("_gp")


def export_vds_to_bgen(vds, interval, output_dir, mean_impute_missing):
    outname = f"gene_{interval.start.contig}_{interval.start.position}_{interval.end.position}"
    print(outname)
    sub_vds = hl.vds.filter_intervals(vds, hl.literal([interval]))
    mt = hl.vds.to_dense_mt(sub_vds)
    mt = mt.annotate_rows(
        rsid=mt.locus.contig
        + ":"
        + hl.str(mt.locus.position)
        + "_"
        + mt.alleles[0]
        + "/"
        + mt.alleles[1]
    )
    mt = mt.annotate_entries(
        GT=hl.if_else(mt.LGT.is_haploid(), hl.call(mt.LGT[0], mt.LGT[0]), mt.LGT)
    )
    mt = gt_to_gp(mt)
    mt = impute_missing_gp(mt, mean_impute=mean_impute_missing)
    hl.export_bgen(mt, f"{output_dir}/{outname}", gp=mt.GP, varid=mt.rsid)


def get_vat_field_types():
    tags = (
        "afr",
        "amr",
        "asj",
        "eas",
        "eur",
        "fin",
        "mid",
        "nfr",
        "sas",
        "oth",
        "max",
        "all",
    )
    fields = ("ac", "an", "sc")
    types = {}
    for tag in tags:
        types[f"gvs_{tag}_af"] = hl.tfloat64
        types[f"gnomad_{tag}_af"] = hl.tfloat64
        for field in fields:
            types[f"gvs_{tag}_{field}"] = hl.tint32
            types[f"gnomad_{tag}_{field}"] = hl.tint32
    types["is_canonical_transcript"] = hl.tbool
    types["omim_phenotypes_id"] = hl.tarray(hl.tint32)
    for x in [
        "position",
        "gene_omim_id",
        "splice_ai_acceptor_gain_distance",
        "splice_ai_acceptor_loss_distance",
        "splice_ai_donor_gain_distance",
        "splice_ai_donor_loss_distance",
    ]:
        types[x] = hl.tint32
    for x in [
        "revel",
        "splice_ai_acceptor_gain_score",
        "splice_ai_acceptor_loss_score",
        "splice_ai_donor_gain_score",
        "splice_ai_donor_loss_score",
    ]:
        types[x] = hl.tfloat64
    for x in [
        "omim_phenotypes_name",
        "clinvar_classification",
        "clinvar_phenotype",
        "consequence",
        "dbsnp_rsid",
    ]:
        types[x] = hl.tarray(hl.tstr)
    return types

