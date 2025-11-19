#!/usr/bin/env python3

__author__ = "wlu"

import argparse

from aou_gwas.utils.utils import *  # for QoB
from aou_gwas.utils.annotations import * # for QoB
from aou_gwas.utils.resources import *  # for QoB
from aou_gwas.utils.resources import annotate_adj  # for QoB
import hail as hl
import hailtop.fs as hfs


def main(args):
    analysis_type = "variant" if args.single_variant_only else "gene"
    mt_name = 'EXOME MT' if analysis_type == 'gene' else 'ACAF MT'
    pops = args.pops.split(",") if (args.pops is not None) else POPS
    app_name = ''
    if args.run_vep:
        app_name = 'run_vep'
    if args.run_relatedness:
        app_name = f'run_relatedness_{args.relatedness_type}'
    if args.get_duplicated_samples:
        app_name = 'run_duplicates'
    if args.run_sample_qc:
        app_name = 'run_sample_qc'
    if args.run_pca:
        app_name = f'run_pca_{args.pops.replace(",","_")}{"_project" if args.project_related_samples else ""}{"_R2" if args.pca_pruned_samples else "_R1"}'
    if args.create_gene_mapping_file:
        app_name = f'create_gene_map_file_{args.pops.replace(",","_")}'
    if args.create_call_stats_ht:
        app_name = f'create_call_stats_ht_{args.pops.replace(",", "_")}{"_pruned" if args.pca_pruned_samples else ""}[{mt_name}]'
    if args.count_bgen_n_var:
        app_name = f'count_bgen_n_var_{args.pops.replace(",", "_")}'
    if args.process_mirna_mt:
        app_name = f'process_miRNA_MT'
    if args.update_saige_intervals:
        app_name = f'update_SAIGE_{analysis_type}_intervals'
    if args.process_brava_annotations:
        app_name = f'Process_BRaVa_annotations'
    hl.stop()
    hl.init(
        tmp_dir='gs://aou_tmp/',
        driver_memory="highmem",
        driver_cores=8,
        worker_memory="highmem",
        worker_cores=1,
        default_reference="GRCh38",
        log=f"~/PycharmProjects/aou_gwas/log/pre_process_saige_data_{app_name}_{date.today()}.log",
        app_name=app_name
    )
    # hl._set_flags(profile='1')

    if args.count_bgen_n_var:
        analysis_type = "variant" if args.single_variant_only else "gene"
        print(f'Analysis type: {analysis_type}')
        def count_var_bgen_from_mt(
                pop,
                analysis_type,
                mean_impute_missing: bool = True,
                no_adj: bool = False,
                variant_ac_filter: int = 0,
                variant_callrate_filter: float = 0,
        ):
            mt = get_filtered_mt(analysis_type=analysis_type, filter_variants=True, filter_samples=True,
                                 adj_filter=not no_adj, pop=pop)

            mt = mt.select_entries("GT")
            mt = mt.filter_rows(
                hl.agg.count_where(mt.GT.is_non_ref()) > 0
            )  # Filter to non-reference sites
            mt = mt.annotate_rows(
                rsid=mt.locus.contig + ":" + hl.str(mt.locus.position) + "_" + mt.alleles[0] + "/" + mt.alleles[1]
            )  # Annotate rsid

            if variant_ac_filter:
                mt = mt.filter_rows(mt.info.AC[0] > variant_ac_filter)

            if variant_callrate_filter:
                mt = mt.filter_rows(mt.info.AN / (2 * N_SAMPLES_RAW[pop]) > variant_callrate_filter)

            mt = mt.annotate_entries(
                GT=hl.if_else(mt.GT.is_haploid(), hl.call(mt.GT[0], mt.GT[0]), mt.GT)
            )
            mt = gt_to_gp(mt)
            mt = impute_missing_gp(mt, mean_impute=mean_impute_missing)
            print(f'TOTAL N VAR in {pop} Bgen: {mt.count()}')

        pops = args.pops.split(",") if (args.pops is not None) else POPS
        for pop in pops:
            count_var_bgen_from_mt(pop=pop, analysis_type = analysis_type)

    if args.run_vep:
        # Vep resource from gnomAD:
        # https://github.com/broadinstitute/gnomad_methods/blob/e2f2b055a2d40991851e7852dea435b4eaef45ea/gnomad/resources/grch38/reference_data.py#L100
        if (not hfs.exists(get_aou_util_path(name="vep"))) or args.overwrite:
            # Load info tables
            # vat_ht = hl.read_table(get_aou_util_path(name="vat"))
            # vat_ht = vat_ht.collect_by_key()
            # vat_ht.write(get_aou_util_path(name="vat_collected_by_key"), overwrite=True)
            # vat_ht = hl.read_table(get_aou_util_path(name="vat_collected_by_key"))
            # indel_ht = vat_ht.filter(~hl.is_snp(vat_ht.alleles[0], vat_ht.alleles[1]))
            # indel_ht.write(get_aou_util_path(name="vat_collected_by_key_indel"), overwrite=True)
            # vat_ht = hl.read_table(get_aou_util_path(name="vat_collected_by_key_indel"))
            # print(vat_ht.count())
            # # #### Code for testing
            # vat_ht = vat_ht.filter((vat_ht.locus.contig == 'chr1') & (vat_ht.locus.position < 10330))
            # # vat_ht = vat_ht.filter((vat_ht.locus.contig == 'chr21'))
            # # vep_ht = hl.vep(vat_ht, csq=False)
            # vep_ht = vep_or_lookup_vep(vat_ht, vep_version="vep105")
            # vep_ht.write('gs://aou_wlu/test/vep_test.ht', overwrite=True)
            # vep_ht = hl.read_table('gs://aou_wlu/test/vep_test.ht')
            # vep_ht.describe()
            # vep_ht.show()
            # print(vep_ht.count())
            #### Run the code below, Jan 11, 2024
            vat_ht = hl.read_table(get_aou_util_path(name="vat_collected_by_key_indel"))
            vep_ht = vep_or_lookup_vep(vat_ht, vep_version="vep105")
            vep_ht.write(
                get_aou_util_path(name="vep"), # vep for indels
                overwrite=args.overwrite,
            )
        vep_ht = hl.read_table(get_aou_util_path(name="vep")) # vep for indels
        vep_ht.describe()
        # vep_ht.show()
        # print(vep_ht.count())
        if (not hfs.exists(f'{get_aou_util_path(name="vep_full")}/_SUCCESS')) or args.overwrite:
            snp_vep_ht = hl.read_table(
                'gs://gcp-public-data--gnomad/resources/context/grch38_context_vep_annotated.v105.ht')
            indel_vep_ht = hl.read_table(get_aou_util_path(name="vep"))

            def process_vep_ht(ht):
                ht = process_consequences(ht)
                ht = ht.explode(ht.vep.worst_csq_by_gene_canonical)
                ht = ht.filter(ht.vep.worst_csq_by_gene_canonical.gene_id.startswith('ENSG'))
                ht = ht.annotate(
                    variant_id=ht.locus.contig + ':' + hl.str(ht.locus.position) + ':' + ht.alleles[0] + ':' +
                               ht.alleles[1],
                    annotation=annotation_case_builder(ht.vep.worst_csq_by_gene_canonical))
                ht = ht.select(
                    annotation = ht.annotation,
                    variant_id = ht.variant_id,
                    gene_id=ht.vep.worst_csq_by_gene_canonical.gene_id,
                    gene_symbol=ht.vep.worst_csq_by_gene_canonical.gene_symbol,
                    transcript_id=ht.vep.worst_csq_by_gene_canonical.transcript_id,
                    polyphen2=ht.vep.worst_csq_by_gene_canonical.polyphen_prediction,
                    amino_acids=ht.vep.worst_csq_by_gene_canonical.amino_acids,
                    lof=ht.vep.worst_csq_by_gene_canonical.lof,
                    most_severe_consequence=ht.vep.worst_csq_by_gene_canonical.most_severe_consequence,
                )
                return ht

            snp_vep_ht = process_vep_ht(snp_vep_ht)
            snp_vep_ht = snp_vep_ht.checkpoint('gs://aou_tmp/vep/snp_vep_250k.ht', _read_if_exists=True)
            indel_vep_ht = process_vep_ht(indel_vep_ht)
            indel_vep_ht = indel_vep_ht.checkpoint('gs://aou_tmp/vep/indel_vep_250k.ht', _read_if_exists=True)
            vep_ht = snp_vep_ht.union(indel_vep_ht)

            vep_ht.checkpoint(get_aou_util_path(name="vep_full"), overwrite = True)

        vep_ht = hl.read_table(get_aou_util_path(name="vep_full"))
        vep_ht.describe()
        vep_ht.show()
        print(vep_ht.count()) # 7037600680

    if args.update_saige_intervals:
        analysis_type = "variant" if args.single_variant_only else "gene"
        print(f'Analysis type: {analysis_type}')

        if analysis_type == 'variant':
            chunk_size = CHUNK_SIZE['afr']
            variant_interval_path = get_saige_interval_path(analysis_type, chunk_size)
            if not hfs.exists(variant_interval_path) or args.overwrite:
                print(f'Generating interval file for SAIGE (chunk size: {chunk_size})...')
                intervals = []
                for chrom in CHROMOSOMES:
                    chromosome = f"chr{chrom}"
                    CHROMOSOME_LENGTHs = hl.get_reference(REFERENCE).lengths
                    chrom_length = CHROMOSOME_LENGTHs[chromosome]
                    for start_pos in range(1, chrom_length, chunk_size):
                        end_pos = (
                            chrom_length
                            if start_pos + chunk_size > chrom_length
                            else (start_pos + chunk_size)
                        )
                        interval = hl.Interval(
                            hl.Locus(chromosome, start_pos, reference_genome='GRCh38'),
                            hl.Locus(chromosome, end_pos, reference_genome='GRCh38')
                        )

                        intervals.append(interval)
                import pandas as pd
                df = pd.DataFrame({'interval': intervals})
                ht = hl.Table.from_pandas(df)
                ht.write(variant_interval_path, overwrite=args.overwrite)
            ht = hl.read_table(variant_interval_path)
            ht.describe()
            ht.show(10)
            print(ht.count())
        else:
            n_genes_per_group = N_GENE_PER_GROUP
            gene_interval_path = get_saige_interval_path(analysis_type, n_genes_per_group)
            if not hfs.exists(gene_interval_path) or args.overwrite:
                print(f'Generating interval file for SAIGE-GENE (N genes per group: {n_genes_per_group})...')
                gtf = hl.experimental.import_gtf(GTF_PATH,
                                                 reference_genome='GRCh38',
                                                 skip_invalid_contigs=True)
                ht = group_gene_interval(gtf, n=n_genes_per_group, path=gene_interval_path, overwrite=args.overwrite)
            ht = hl.read_table(gene_interval_path)
            ht.describe()
            ht = ht.annotate(n_genes = hl.len(ht.genes))
            ht.filter(ht.n_genes > N_GENE_PER_GROUP).show(100)
            ht.filter(ht.n_genes < N_GENE_PER_GROUP/2).show(100)
            print(ht.aggregate(hl.agg.counter(ht.n_genes)))
            print(ht.count())

    if args.create_call_stats_ht:
        pruned = args.pca_pruned_samples
        pops.append('all')
        for pop in pops:
            path = get_call_stats_ht_path(pop=pop, pruned=pruned, analysis_type=analysis_type)
            if not hfs.exists(path):
                print(f'[{mt_name}] Exporting sample ({"pruned" if pruned else "raw"}) call stats ({pop.upper()})')
                print(f'[{mt_name}] To path: {path}')
                mt = get_filtered_mt(analysis_type=analysis_type, filter_variants=True, filter_samples=True,
                                     adj_filter=True, pop=pop, prune_samples=pruned)
                call_stats_ht = mt.annotate_rows(call_stats=hl.agg.call_stats(mt.GT, mt.alleles)).rows()
                call_stats_ht = call_stats_ht.naive_coalesce(1000).checkpoint(path)
        path = get_call_stats_ht_path(pop='full', pruned=pruned, analysis_type=analysis_type)
        if not hfs.exists(f'{path}/_SUCCESS') or args.overwrite:
            print(f'[{mt_name}] Exporting FULL call stats to path: {path}')
            ht = hl.read_table(get_aou_util_path(name="vep_full"))
            global_ht = hl.read_table(get_call_stats_ht_path(pop='all', pruned=False, analysis_type=analysis_type))
            ht = ht.filter(hl.is_defined(global_ht[ht.key]))
            ht = ht.annotate(
                freq = hl.struct(**{f'{pop.upper()}': hl.read_table(get_call_stats_ht_path(pop=pop, pruned=True, analysis_type=analysis_type))[ht.key].call_stats for pop in pops}),
                freq_raw = hl.struct(**{f'{pop.upper()}': hl.read_table(get_call_stats_ht_path(pop=pop, pruned=False, analysis_type=analysis_type))[ht.key].call_stats for pop in pops})
            )
            ht.describe()
            ht.naive_coalesce(1000).checkpoint(path, overwrite=args.overwrite)
        ht = hl.read_table(path)
        ht.describe()
        ht.show()
        print(ht.count())

    if args.compute_caf:
        ht = hl.read_table(get_call_stats_ht_path(pop='all', pruned=True, analysis_type='gene'))
        ht = ht.filter(hl.is_defined(ht.freq))
        ht = ht.annotate(annotation = hl.if_else(hl.literal({"missense", "LC"}).contains(ht.annotation), "missenseLC", ht.annotation))
        sub_ht = ht.filter(
            hl.literal({"missenseLC", "pLoF"}).contains(ht.annotation)
        )
        sub_ht = sub_ht.key_by()
        sub_ht = sub_ht.annotate(annotation="pLoF;missenseLC")
        for pop in POPS:
            print(pop.upper())
            caf_ht = ht.group_by('gene_id', 'gene_symbol', 'annotation').aggregate(**{ f'CAF_{pop}':hl.agg.sum(ht.freq[pop.upper()].AF[1])},
                                                                                   **{f'CAF_{pop}_0.001': hl.agg.filter(ht.freq[pop.upper()].AF[1] < 0.001 ,hl.agg.sum(ht.freq[pop.upper()].AF[1]))})
            sub_caf_ht = sub_ht.group_by('gene_id', 'gene_symbol', 'annotation').aggregate(**{f'CAF_{pop}':hl.agg.sum(sub_ht.freq[pop.upper()].AF[1])},
                                                                                           **{f'CAF_{pop}_0.001': hl.agg.filter(sub_ht.freq[pop.upper()].AF[1] < 0.001 ,hl.agg.sum(sub_ht.freq[pop.upper()].AF[1]))})
            caf_ht = caf_ht.union(sub_caf_ht)
            caf_ht.describe()
            caf_ht.write(get_call_stats_ht_path(pop=pop, pruned=True, analysis_type='caf'), overwrite=True)

    if args.create_gene_mapping_file:
        snp_vep_ht = hl.read_table('gs://gcp-public-data--gnomad/resources/context/grch38_context_vep_annotated.v105.ht')
        snp_vep_ht = snp_vep_ht.drop('grange', 'vep_help', 'vep_config', 'version')
        indel_vep_ht = hl.read_table(get_aou_util_path(name="vep"))
        for pop in pops:
            if not args.skip_raw_gene_map_file:
                call_stats_ht = hl.read_table(get_call_stats_ht_path(pop=pop, pruned=True, analysis_type='gene'))
                call_stats_ht = call_stats_ht.filter(call_stats_ht.call_stats.AC[1]>0)
                print(f'---------Generating raw gene mapping HT ({pop.upper()})-----------------')
                max_an = call_stats_ht.aggregate(
                    hl.struct(
                        autosomes=hl.agg.max(call_stats_ht.call_stats.AN),
                        x=hl.agg.filter(
                            call_stats_ht.locus.in_x_nonpar(),
                            hl.agg.max(call_stats_ht.call_stats.AN)
                        ),
                        y=hl.agg.filter(
                            call_stats_ht.locus.in_y_nonpar(),
                            hl.agg.max(call_stats_ht.call_stats.AN),
                        ),
                    ),
                )

                an = call_stats_ht.call_stats.AN
                call_stats_ht = call_stats_ht.filter(
                    hl.case()
                    .when(call_stats_ht.locus.in_x_nonpar(), an > 0.8 * max_an.x)
                    .when(call_stats_ht.locus.in_y_nonpar(), an > 0.8 * max_an.y)
                    .default(an > 0.8 * max_an.autosomes)
                )

                snp_vep_ht = snp_vep_ht.annotate(freq= call_stats_ht[snp_vep_ht.key].call_stats.AF[1])
                snp_vep_ht = snp_vep_ht.filter(
                    hl.is_defined(snp_vep_ht.freq)
                )
                indel_vep_ht = indel_vep_ht.annotate(freq=call_stats_ht[indel_vep_ht.key].call_stats.AF[1])
                indel_vep_ht = indel_vep_ht.filter(
                    hl.is_defined(indel_vep_ht.freq)
                )
                gene_map_ht = create_gene_map_ht(snp_vep_ht, indel_vep_ht, freq_field='freq')
                print(f'---------Exporting raw gene mapping HT ({pop.upper()})-----------------')
                gene_map_ht.checkpoint(
                    get_aou_gene_map_ht_path(pop=pop, processed=False), _read_if_exists=True
                )

            gene_map_ht = hl.read_table(
                get_aou_gene_map_ht_path(pop=pop, processed=False),
            )
            gene_map_ht.describe()

            gene_map_ht = post_process_gene_map_ht(gene_map_ht, freq_cutoff=0.01)
            print(f'---------Exporting processed gene mapping HT ({pop.upper()})-----------------')
            print(get_aou_gene_map_ht_path(pop=pop, processed=True))
            gene_map_ht = gene_map_ht.checkpoint(
                get_aou_gene_map_ht_path(pop=pop, processed=True), _read_if_exists=not args.overwrite, overwrite = args.overwrite
            )
            gene_map_ht.describe()
            gene_map_ht.show()

    if args.get_duplicated_samples:
        # TODO: follow up on zulip https://hail.zulipchat.com/#narrow/stream/123010-Hail-Query-0.2E2-support/topic/ClassCastException/near/394696118
        if (not hfs.exists(get_aou_relatedness_path(extension="duplicates.ht"))) or args.overwrite:
            print('Generating duplicated samples...')
            degree_0th_cutoff = 0.354
            ht = hl.read_table(get_aou_util_path(name="relatedness"))
            ht = ht.filter(ht.kin > degree_0th_cutoff)
            print('--------------------Running hl.maximal_independent_set()----------------')
            duplicated_samples_to_remove = hl.maximal_independent_set(ht['i.s'], ht['j.s'], keep=False)
            duplicated_samples_to_remove.describe()
            duplicated_samples_to_remove.checkpoint(get_aou_relatedness_path(extension="duplicates.ht"),
                                                    overwrite=args.overwrite,
                                                    _read_if_exists= not args.overwrite)
        duplicated_samples_to_remove = hl.read_table(get_aou_relatedness_path(extension="duplicates.ht"))
        print(f'N duplicated samples: {duplicated_samples_to_remove.count()}')

        if (not hfs.exists(get_aou_relatedness_path(extension="1st_degrees.ht"))) or args.overwrite:
            print('Generating 0th degree samples...')
            ht = hl.read_table(get_aou_util_path(name="relatedness", parsed=True))
            ht = ht.filter(ht.kin > 0.354)
            ht = ht.key_by()
            ht_i = ht.key_by(s=ht['i.s']).select()
            ht_j = ht.key_by(s=ht['j.s']).select()
            samples_1st_degree = ht_i.union(ht_j).distinct()
            samples_1st_degree.describe()
            samples_1st_degree.checkpoint(get_aou_relatedness_path(extension="1st_degrees.ht"),
                                                    overwrite=args.overwrite,
                                                    _read_if_exists= not args.overwrite)
        samples_1st_degree = hl.read_table(get_aou_relatedness_path(extension="1st_degrees.ht"))
        print(f'N duplicated samples: {samples_1st_degree.count()}')

    if args.run_relatedness:
        MAF_CUTOFF = 0.05
        mt = hl.read_matrix_table(ACAF_MT_PATH)  # (99250816, 245394)
        mt.describe()
        # print(mt.info.AF.summarize())
        # # - AF (array<float64>):
        # #   Non-missing: 99250796 (100.00%)
        # #       Missing: 20 (0.00%)
        # #      Min Size: 1
        # #      Max Size: 1
        # #     Mean Size: 1.00
        # #
        # #   - AF[<elements>] (float64):
        # #     Non-missing: 99250796 (100.00%)
        # #         Missing: 0
        # #         Minimum: 0.00
        # #         Maximum: 1.00
        # #            Mean: 0.04
        # #         Std Dev: 0.13
        # print(mt.aggregate_rows(hl.agg.count_where(mt.info.AF[0] > MAF_CUTOFF))) # 10289329
        ### Filter to 1M random common variants before running IBD

        if (
                not hfs.exists(get_aou_relatedness_path(extension="100Kvar.ht"))
        ) or args.overwrite:
            print(f"-------------Downsampling to 100K common variants-----------")
            ht = hl.read_table('gs://aou_analysis/250k/data/aou_ACAF_rows_full.ht')
            sampled_variants = ht.aggregate(
                hl.agg.filter(
                    ht.info.AF[0] >= MAF_CUTOFF,
                    hl.agg._reservoir_sample(ht.key, 100000),
                )
            )
            variants = [variant for variant in sampled_variants]
            ht = hl.Table.parallelize(variants).key_by(*ht.key.keys())
            ht = ht.checkpoint(
                get_aou_relatedness_path(extension="100Kvar.ht"),
                _read_if_exists=(not args.overwrite),
                overwrite=args.overwrite,
            )
        ht = hl.read_table(get_aou_relatedness_path(extension="100Kvar.ht"))
        # print(f'N variants: {ht.count()}')
        mt = mt.filter_rows(hl.is_defined(ht[mt.row_key]))
        mt = mt.naive_coalesce(1000).checkpoint(get_aou_relatedness_path(extension="100Kvar.mt"))

        if args.relatedness_type == "ibd":
            print(f"-------------Running hl.identity_by_descent() -----------")
            mt = mt.unfilter_entries()
            mt = mt.annotate_entries(mean_gt = mt.aggregate_entries(hl.agg.mean(mt.GT.n_alt_alleles(), _localize=False)))
            mt = mt.annotate_entries(GT = hl.or_missing(mt.GT, mt.mean_gt))
            relatedness_ht = hl.identity_by_descent(mt, maf=mt.info.AF[0])
            if args.overwrite or (
                    not hfs.exists(f'{get_aou_relatedness_path(extension="ht")}')
            ):
                print(f"-------------Writing AoU IBD HT -----------")
                relatedness_ht.write(
                    get_aou_relatedness_path(extension="ibd.ht"), args.overwrite
                )
        if args.relatedness_type == "king":
            print(f"-------------Running hl.king() -----------")
            relatedness_mt = hl.king(mt.GT)
            if args.overwrite or (
                    not hfs.exists(f'{get_aou_relatedness_path(extension="king.ht")}')
            ):
                print(f"-------------Writing AoU King relatedness MT -----------")
                relatedness_mt.write(
                    get_aou_relatedness_path(extension="king.mt"), args.overwrite
                )
    if args.run_sample_qc:
        if not hfs.exists(get_aou_util_path('mt_sample_qc')):
            print('Run sample qc MT.....')
            mt = hl.read_matrix_table(ACAF_MT_PATH)
            mt = mt.filter_rows(mt.locus.in_autosome())
            # mt = mt.filter_rows(mt.locus.contig == 'chr1')
            ht = hl.sample_qc(mt, name='mt_sample_qc')
            ht.write(get_aou_util_path('mt_sample_qc'), overwrite=args.overwrite)
        mt = hl.read_matrix_table(get_aou_util_path('mt_sample_qc'))
        ht = mt.cols()
        print('Write sample QC HT.....')
        ht.write(get_aou_util_path('acaf_mt_sample_qc'))

    if args.run_pca:
        meta_ht = hl.read_table(get_sample_meta_path(annotation=True))
        meta_ht = meta_ht.filter(~meta_ht.related)
        pops = args.pops.split(",") if (args.pops is not None) else POPS
        tag = "pruned_" if args.pca_pruned_samples else ""
        for pop in pops:
            if not args.dataproc and not args.project_related_samples:
                call_stats_ht = hl.read_table(get_call_stats_ht_path(pop=f'{pop}_repartition' if pop in ['afr', 'amr'] else pop, pruned=False, analysis_type='variant'))
                ### Created from code below:
                # pops = ['afr', 'mid', 'sas', 'amr', 'eas', 'eur']
                # for pop in pops:
                #     print(f'--------------{pop}--------------')
                #     mt = get_filtered_mt(analysis_type='variant', filter_variants=True, filter_samples=True,
                #                          adj_filter=True, pop=pop, prune_samples=False)
                #     call_stats_ht = mt.annotate_rows(call_stats=hl.agg.call_stats(mt.GT, mt.alleles)).rows()
                #     call_stats_ht = call_stats_ht.select('call_stats')
                #     call_stats_ht.naive_coalesce(1000).checkpoint(get_call_stats_ht_path(pop=pop, pruned=False), _read_if_exists=True)
                #     f"{DATA_PATH}/utils/call_stats/aou_{pop}_variant_info_{tranche}.{extension}"
                #     call_stats_ht.describe()
                ############################
                mt = get_filtered_mt(analysis_type='variant', filter_variants=True, filter_samples=True,
                                     adj_filter=True, pop=pop, prune_samples=False)
                print("------Removing HLA & inversion regions------")
                # Common inversion taken from Table S4 of https://www.ncbi.nlm.nih.gov/pubmed/27472961
                # (converted to GRCh38 by: https://liftover.broadinstitute.org/#input=chr8%3A8055789-11980649&hg=hg19-to-hg38 )
                # Also removing HLA, from https://www.ncbi.nlm.nih.gov/grc/human/regions/MHC?asm=GRCh38
                mt = mt.filter_rows(
                    ~hl.parse_locus_interval(
                        "chr8:8198267-12123140", reference_genome="GRCh38"
                    ).contains(mt.locus)
                    & ~hl.parse_locus_interval(
                        "chr6:28510120-33480577", reference_genome="GRCh38"
                    ).contains(mt.locus)
                )
                MIN_CR = CALLRATE_CUTOFF
                print(f'------Call rate filter: {MIN_CR}------')
                MIN_AF = 0.01
                print(f'------AF filter: {MIN_AF}------')
                print(f'------Number of samples: {N_SAMPLES_HARD_FILTERED[pop]}------')
                variants_to_keep = call_stats_ht.filter(
                    (call_stats_ht.locus.in_autosome()) &
                    (hl.is_snp(call_stats_ht.alleles[0], call_stats_ht.alleles[1])) &
                    (call_stats_ht.call_stats.AF[1] >= MIN_AF) &
                    ((call_stats_ht.call_stats.AN >= (N_SAMPLES_HARD_FILTERED[pop] * 2 * MIN_CR)))
                )
                print('Filtering Variants...')
                mt = mt.filter_rows(hl.is_defined(variants_to_keep[mt.row_key])) # filter to high quality variants
                print('Filtering Samples...')
                mt = mt.filter_cols(hl.is_defined(meta_ht[mt.col_key])) # filter to unrelated samples -> later to project
                mt = mt.unfilter_entries()
                mt = mt.naive_coalesce(250).checkpoint(f'gs://aou_analysis/250k/data/utils/pca/pca_{pop}.mt',
                                                       _read_if_exists=not args.overwrite, overwrite=args.overwrite)

                mt = hl.read_matrix_table(f'gs://aou_analysis/250k/data/utils/pca/pca_{pop}.mt')
                print(f'-----------{pop.upper()}-----------')
                print(mt.count())

            if args.dataproc:
                if args.pca_mt_prep:
                    mt = hl.read_matrix_table(f'gs://aou_analysis/250k/data/utils/pca/pca_{pop}.mt')
                    print(f'-----------{pop.upper()}-----------')
                    print(mt.count())
                    # Downsample to 1M variants -> LD pruning -> filter MT to LD pruned variants -> repartition MT to 50
                    prop = round(1000000 / mt.count_rows(), 5)
                    print(f'Proportion of variants: {prop}')
                    print(f'Number of variants: {prop * mt.count_rows()}')
                    mt = mt.sample_rows(prop)
                    mt = mt.checkpoint(f'gs://aou_analysis/250k/data/utils/pca/pca_{pop}_downsampled.mt',
                                                               _read_if_exists=not args.overwrite, overwrite=args.overwrite)

                    print(f"-------------Exporting the LD pruned downsampled variant HT (pop: {pop})-----------")
                    # TODO:
                    # 1) run ld_prune
                    # 2) HGDP + 1KG ld pruned variants
                    ht = hl.ld_prune(mt.GT,
                                     r2=0.1,
                                     bp_window_size=int(1e7),
                                     block_size=1024,
                                     )
                    ht = ht.checkpoint(
                        f'gs://aou_analysis/250k/data/utils/pca/pca_{pop}_ld_pruned.ht',
                        _read_if_exists=not args.overwrite,
                        overwrite=args.overwrite,
                    )

                    mt = mt.filter_rows(hl.is_defined(ht[mt.row_key]))
                    mt = mt.checkpoint(f'gs://aou_analysis/250k/data/utils/pca/pca_{pop}_ld_pruned.mt',
                                                          _read_if_exists=not args.overwrite, overwrite=args.overwrite)



                    mt.naive_coalesce(50).checkpoint(f'gs://aou_analysis/250k/data/utils/pca/pca_{pop}_50partitions.mt',
                                                          _read_if_exists=not args.overwrite, overwrite=args.overwrite)

                print(f'-----------{pop.upper()} LD-pruned & repartitioned MT-----------')
                mt = hl.read_matrix_table(f'gs://aou_analysis/250k/data/utils/pca/pca_{pop}_50partitions.mt')
                print(mt.count())  # per population
                if args.pca_pruned_samples:
                    print('Pruning samples...')
                    pruned_ht = hl.import_table(CENTROID_PRUNED_SAMPLES, types = {'s':hl.tstr}).key_by('s')
                    mt = mt.filter_cols(hl.is_defined(pruned_ht[mt.col_key]))
                    mt = mt.checkpoint(f'gs://aou_analysis/250k/data/utils/pca/pca_{pop}_sample_pruned.mt',
                                       _read_if_exists=not args.overwrite, overwrite=args.overwrite)
                    print(mt.count())

                print('Running PCA...')
                eigenvalues, scores, loadings = hl.hwe_normalized_pca(
                    mt.GT,
                    compute_loadings=True,
                    k=50,
                )

                print('Writing tables...')
                scores.write(
                    get_pca_ht_path(pop=pop, name=f'{tag}scores'),
                    overwrite=args.overwrite,
                )
                loadings.write(
                    get_pca_ht_path(pop=pop, name=f'{tag}loadings'),
                    overwrite=args.overwrite,
                )

                scores_ht = hl.read_table(get_pca_ht_path(pop=pop, name=f'{tag}scores'))
                scores_ht.describe()

                loadings_ht = hl.read_table(get_pca_ht_path(pop=pop, name=f'{tag}loadings'))
                loadings_ht.describe()

            if args.project_related_samples:
                if args.pca_pruned_samples:
                    print('Use pruned samples...')
                    pca_mt = hl.read_matrix_table(f'gs://aou_analysis/250k/data/utils/pca/pca_{pop}_sample_pruned.mt')
                    mt = get_filtered_mt(analysis_type='variant', filter_variants=True, filter_samples=True,
                                         adj_filter=True, pop=pop, prune_samples=True)
                    # print(mt.count())
                else:
                    print('Use raw samples...')
                    pca_mt = hl.read_matrix_table(f'gs://aou_analysis/250k/data/utils/pca/pca_{pop}_50partitions.mt')
                    mt = get_filtered_mt(analysis_type='variant', filter_variants=True, filter_samples=True,
                                         adj_filter=True, pop=pop, prune_samples=False)
                pca_mt = pca_mt.annotate_rows(pca_af=hl.agg.mean(pca_mt.GT.n_alt_alleles()) / 2)

                mt = mt.filter_cols(hl.is_missing(meta_ht[mt.col_key]))

                pca_loadings = hl.read_table(get_pca_ht_path(pop=pop, name=f'{tag}loadings'))
                pca_loadings = pca_loadings.annotate(
                    pca_af=pca_mt.rows()[pca_loadings.key].pca_af
                )
                pca_ht = hl.experimental.pc_project(
                    mt.GT,
                    pca_loadings.loadings,
                    pca_loadings.pca_af,
                )

                pca_ht.checkpoint(
                    get_pca_ht_path(pop=pop, name=f'{tag}related_scores'),
                    overwrite=args.overwrite,
                )

                related_scores_ht = hl.read_table(get_pca_ht_path(pop=pop, name=f'{tag}related_scores'))
                print(f"Number of related samples : {related_scores_ht.count()}")
                scores_ht = hl.read_table(get_pca_ht_path(pop=pop, name=f'{tag}scores'))
                full_scores_ht = related_scores_ht.union(scores_ht)
                full_scores_ht = full_scores_ht.annotate(**{f"PC{i + 1}": full_scores_ht.scores[i] for i in range(50)})

                full_scores_ht.checkpoint(
                    get_pca_ht_path(pop=pop, name=f'{tag}full_scores'),
                    overwrite=args.overwrite,
                )
                full_scores_ht.describe()

                print(full_scores_ht.count())

                full_scores_ht.export(get_pca_ht_path(pop=pop, name=f'{tag}full_scores', extension='txt.bgz'))

    if args.process_mirna_mt:
        miRNA_ht = hl.read_table(get_aou_util_path(name='miRNA_info'))
        miRNA_ht.describe()
        miRNA_ht.show()
        # mt = hl.read_matrix_table(get_aou_util_path(name='miRNA_raw', extension='mt'))
        # mt.describe()
        # print(mt.count())
        # print('----------Multi-splitting miRNA MT-------------')
        # mt = hl.split_multi(mt)
        # mt = mt.checkpoint(get_aou_util_path(name='miRNA_split', extension='mt'), _read_if_exists = True)
        # mt = mt.naive_coalesce(10).checkpoint(get_aou_util_path(name='miRNA_split_repartitioned', extension='mt'), _read_if_exists = True)

        print('----------Filtering multi-split miRNA MT-------------')
        mt= hl.read_matrix_table(get_aou_util_path(name='miRNA_split_repartitioned', extension='mt'))
        print(mt.count())
        # print(mt.aggregate_rows(hl.agg.counter(hl.len(mt.alleles))))
        print(mt.aggregate_entries(hl.agg.filter(mt.GT.contains_allele(3), hl.agg.take(hl.tuple([mt.s, mt.locus]), 5))))
        mt = mt.annotate_rows(call_stats=hl.agg.call_stats(mt.GT, mt.alleles))

        # mt = mt.annotate_rows(miRNA=miRNA_ht.index(mt.locus, all_matches=True).id)
        mt = mt.filter_entries(mt.FT | hl.is_missing(mt.FT))
        mt = mt.annotate_rows(call_stats=hl.agg.call_stats(mt.GT, mt.alleles))
        print(mt.aggregate_rows(hl.agg.counter(hl.len(mt.call_stats.AC))))
        # mt = mt.annotate_rows(filter_ac = ((mt.call_stats.AC[1] > 0) & (hl.len(mt.call_stats.AC)>1)) | (hl.len(mt.call_stats.AC)<1))
        # mt = mt.filter_rows(mt.filter_ac)
        mt = mt.filter_rows(((mt.call_stats.AC[1] > 0) & (hl.len(mt.call_stats.AC)==2)))
        mt = annotate_adj(mt)
        mt = mt.filter_entries(mt.adj)
        mt = mt.filter_rows(hl.agg.any(mt.GT.n_alt_alleles() >0))
        mt.describe()
        mt = mt.checkpoint(get_aou_util_path(name='miRNA_split_filtered', extension='mt'), overwrite=True)
        mt = mt.annotate_rows(miRNA=miRNA_ht.index(mt.locus, all_matches=True).id)
        mt.describe()

        ht = mt.rows()
        ht.describe()
        ht.show()
        ht = ht.explode(ht.miRNA)
        ht = ht.checkpoint(get_aou_util_path(name='miRNA_split', extension='ht'), overwrite=True)
        print(ht.aggregate(hl.agg.counter(ht.miRNA)))
        print(ht.aggregate(hl.agg.filter(ht.call_stats.AC[1] > 0, hl.agg.counter(ht.miRNA))))

    if args.process_brava_annotations:
        vat_ht = hl.read_table(get_aou_util_path(name='vat'))
        vat_ht.describe()
        print(vat_ht.count())

        def clean_vep_ht(vep_ht):
            from gnomad.utils.vep import process_consequences
            vep_ht = process_consequences(vep_ht)
            vep_ht = vep_ht.explode(vep_ht.vep.worst_csq_by_gene_canonical)
            vep_ht = vep_ht.select(
                lof=vep_ht.vep.worst_csq_by_gene_canonical.lof,
                hgvsp=vep_ht.vep.worst_csq_by_gene_canonical.hgvsp,
                most_severe_csq_variant=vep_ht.vep.worst_csq_for_variant_canonical.most_severe_consequence,
                most_severe_csq_gene=vep_ht.vep.worst_csq_by_gene_canonical.most_severe_consequence,
            )
            return vep_ht

        indel_vep_ht = hl.read_table(get_aou_util_path(name="vep"))
        snp_vep_ht = hl.read_table(
            'gs://gcp-public-data--gnomad/resources/context/grch38_context_vep_annotated.v105.ht')
        indel_vep_ht = clean_vep_ht(indel_vep_ht)
        snp_vep_ht = clean_vep_ht(snp_vep_ht)
        vep_ht = snp_vep_ht.union(indel_vep_ht).distinct()
        vep_ht = vep_ht.checkpoint('gs://aou_wlu/250k_brava/vep_cleaned.ht', _read_if_exists=True)

        def clean_gnomad_ht(gnomad_ht):
            gnomad_ht = gnomad_ht.select(
                cadd=gnomad_ht.in_silico_predictors.cadd,
                gnomad_revel=gnomad_ht.in_silico_predictors.revel_max,
                gnomad_spliceai=gnomad_ht.in_silico_predictors.spliceai_ds_max
            )
            return gnomad_ht

        v4_genome_ht = hl.read_table('gs://gcp-public-data--gnomad/release/4.0/ht/genomes/gnomad.genomes.v4.0.sites.ht')
        v4_exome_ht = hl.read_table('gs://gcp-public-data--gnomad/release/4.0/ht/exomes/gnomad.exomes.v4.0.sites.ht')
        v4_genome_ht = clean_gnomad_ht(v4_genome_ht)
        v4_exome_ht = clean_gnomad_ht(v4_exome_ht)
        v4_ht = v4_genome_ht.union(v4_exome_ht).distinct()
        v4_ht = v4_ht.checkpoint('gs://aou_wlu/250k_brava/gnomad_v4_cleaned.ht', _read_if_exists=True)

        sub_vat = vat_ht.select(
            revel=vat_ht.revel,
            **{f'{field}': vat_ht[field] for field in list(vat_ht.row_value) if field.startswith('splice_ai_')},
            **vep_ht[vat_ht.key],
            **v4_ht[vat_ht.key]
        )

        sub_vat = sub_vat.checkpoint('gs://aou_wlu/250k_brava/variant_annotation_cleaned.ht', _read_if_exists=True)
        sub_vat.describe()

        MISSENSE_CSQS = ["stop_lost", "start_lost", "inframe_insertion", "inframe_deletion", "missense_variant"]
        ht = hl.read_table('gs://aou_wlu/250k_brava/variant_annotation_cleaned.ht')
        ht.describe()
        print(ht.count())
        ht = ht.annotate(splice_ai_ds=hl.max(
            ht.splice_ai_acceptor_gain_score,
            ht.splice_ai_acceptor_loss_score,
            ht.splice_ai_donor_gain_score,
            ht.splice_ai_donor_loss_score,
        ))
        ht = ht.annotate(brava_mask=
                         hl.case(missing_false=True)
                         .when(ht.lof == 'HC', 'High confidence pLoF')
                         .when((hl.literal(MISSENSE_CSQS).contains(ht.most_severe_csq_variant) & (
                                     (ht.revel >= 0.773) | (ht.cadd.phred >= 28.1))) |
                               (ht.splice_ai_ds >= 0.2) |
                               (ht.lof == 'LC'), 'Damaging missense/protein-altering')
                         .when(hl.literal(MISSENSE_CSQS).contains(ht.most_severe_csq_variant),
                               'Other missense/protein-altering')
                         .when((ht.most_severe_csq_variant == 'synonymous_variant') & (ht.splice_ai_ds < 0.2),
                               'Synonymous')
                         .or_missing())
        ht = ht.checkpoint('gs://aou_wlu/250k_brava/variant_annotation_cleaned_brava_annotated.ht', _read_if_exists=True)
        print(ht.aggregate(hl.agg.counter(ht.brava_mask)))



if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--single-variant-only", help="Run single variant SAIGE", action="store_true"
    )
    parser.add_argument(
        "--test",
        help="Whether to run test on a subset of the data (e.g. chr1)",
        action="store_true",
    )
    parser.add_argument(
        "--dataproc",
        help="Whether to run the pipeline on dataproc",
        action="store_true",
    )
    parser.add_argument("--overwrite", help="Overwrite everything", action="store_true")
    parser.add_argument(
        "--pops",
        help="Comma-separated list of pops to run",
        default="afr,amr,eas,eur,mid,sas",
    )
    parser.add_argument(
        "--count-bgen-n-var", help="Count number of variants that gets into bgen file", action="store_true"
    )
    parser.add_argument(
        "--create-call-stats-ht", help="Whether to call stats HT", action="store_true"
    )
    parser.add_argument(
        "--create-gene-mapping-file", help="Whether to create gene mapping file", action="store_true"
    )
    parser.add_argument(
        "--process-brava-annotations", help="Pull annotations for BRaVa masking", action="store_true"
    )
    parser.add_argument(
        "--skip-raw-gene-map-file", help="Whether to skip creating raw gene mapping file", action="store_true"
    )
    parser.add_argument(
        "--run-relatedness",
        help="Compute relatedness information",
        action="store_true",
    )
    parser.add_argument(
        "--get-duplicated-samples",
        help="Get duplicated samples from the original relatedness kinship table",
        action="store_true",
    )
    parser.add_argument(
        "--relatedness-type",
        help="What function to use for running relatedness",
        choices=["ibd", "king"],
    )
    parser.add_argument(
        "--run-sample-qc",
        help="Whether to run sample QC on the ACAF MT",
        action="store_true",
    )
    parser.add_argument(
        "--run-pca",
        help="Whether to run pca",
        action="store_true",
    )
    parser.add_argument(
        "--project-related-samples",
        help="Whether to run PCA projection for the related samples",
        action="store_true",
    )
    parser.add_argument(
        "--pca-pruned-samples",
        help="Whether to run PCA on the pruned samples",
        action="store_true",
    )
    parser.add_argument(
        "--pca-mt-prep",
        help="Whether to re generate the PCA MT",
        action="store_true",
    )
    parser.add_argument(
        "--process-mirna-mt",
        help="Whether to process the miRNA MT",
        action="store_true",
    )
    parser.add_argument(
        "--update-saige-intervals",
        help="Whether to update SAIGE intervals",
        action="store_true",
    )
    parser.add_argument(
        "--compute-caf",
        help="Whether to compute CAF",
        action="store_true",
    )
    parser.add_argument("--run-vep", help="Run Vep on the VAT", action="store_true")
    args = parser.parse_args()
    main(args)
