import argparse
import logging
import hailtop.batch as hb
import pandas as pd
import hailtop.fs as hfs
# from aou_gwas import *  # for dataproc
from utils.utils import *  # for QoB
from utils.resources import *  # for QoB
from utils.results_loading import * # for QoB
from tqdm import tqdm


############################ functions copy-pasted from else where, double check for any updates ############################

TRANCHE = "v8"
MY_BUCKET = 'gs://aou_wlu'
ANALYSIS_BUCKET = "gs://aou_analysis/v8"
DATA_PATH = f'{ANALYSIS_BUCKET}/data'
TMP_BUCKET = 'gs://aou_tmp/v8'


PILOT_PHENOTYPES = ['heart-rate-mean', '3004501', '3022192', 'DE_688.1', 'EM_247']
PILOT_PHENOTYPES = ['EM_236.1', 'CV_401.1', 'A10BJ']
# PILOT_PHENOTYPES = ['EM_202.2', '3010457', 'DE_688.1', 'EM_247', 'weight', 'CV_401', 'EM_239', '3004501', 'blood-pressure-systolic-mean', "height", "heart-rate-mean", "A10BJ", "A10BJ06", "random_0.5_continuous_1", "random_0.5_0.01_1", "random_0.5_0.5_1", "random_0.5_0.2_1", 
#                      "random_0.5_0.1_1", "random_0.5_0.001_1", "3004410_median", "3028288_median", "411", "250.2", "CV_404.1", "836786"]
# PILOT_PHENOTYPES = ['EM_202.2', 'random_0.5_continuous_3', 'random_0.5_0.01_3', '3022192']
# PILOT_PHENOTYPES = ['204.1', '189.21', '195.1', '335', 'NS_326.1', '450', '384.4', 'SO_391.7', '425.1', '260', 
#                      'EM_230.1', '361', '244.1', 'EM_200.11', 'CV_416.2', '174.1', 'RE_469', 'RE_469.2', 'CM_761.13', 'BI_170', 'EM_202', 'GI_513.2', 'EM_256.72'] + ['427.2', '174', '497', '496.2', '427.4', '585.3', '751.21', '697', '288.1', '250', '703', '531.2', '276.5', '275', '560', '280']
# PILOT_PHENOTYPES = ['3010457', 'DE_688.1', 'EM_247', 'weight', 'CV_401', 'EM_239', '3004501', 'blood-pressure-systolic-mean', "height", "heart-rate-mean", "A10BJ", "A10BJ06", "3022192_median", "random_0.5_continuous_1", "random_0.5_0.01_1", "random_0.5_0.5_1", "random_0.5_0.2_1", 
#                      "random_0.5_0.1_1", "random_0.5_0.001_1", "3004410_median", "3028288_median", "411", "250.2", "CV_404.1", "836786"]
# PILOT_PHENOTYPES = ['3010457', 'DE_688.1', 'EM_247', 'weight']
# PILOT_PHENOTYPES = ['3010457_male', '3010457_female', 'DE_688.1_male', 'DE_688.1_female', 'EM_247_male', 'EM_247_female', 'weight_male', 'weight_female']
# PILOT_PHENOTYPES = ["height", "heart-rate-mean", "A10BJ", "A10BJ06", "3022192_median", "random_0.5_continuous_1", "random_0.5_0.01_1", "random_0.5_0.5_1", "random_0.5_0.2_1", 
#                     "random_0.5_0.1_1", "random_0.5_0.001_1", "3004410_median", "3028288_median", "411", "250.2", "CV_404.1", "836786"] 
# PILOT_PHENOTYPES = PILOT_PHENOTYPES + [f'{pheno}_male' for pheno in PILOT_PHENOTYPES] + [f'{pheno}_female' for pheno in PILOT_PHENOTYPES] + \
#                     ["Birth_to_PULMHEART_amr_ALL", "Birth_to_CAD_eas_ALL", "Birth_to_CAD_sas_ALL", "Birth_to_MI_eur_ALL", "Birth_to_DEMENTIA_eur_ALL",
#                     "Birth_to_ATHSCLE_eur_ALL", "HYPTENSESS_to_MI_eur_ALL", "PARKINSON_to_DEMENTIA_eur_ALL", "T2D_to_ATHSCLE_eur_ALL"]


def read_pickle_dict(dict_path: str):
    dict = pd.read_pickle(dict_path)
    return dict

def write_pickle_dict(output: str, dict: dict):
    with hfs.open(output, "wb") as f:
        pickle.dump(dict, f)
    f.close()


#############################################################################################################################

def export_multi_pop_txt_results(analysis_type:str, phenoname:str):
    all_path = f'{get_aou_saige_results_root(analysis_type=analysis_type, pop="all", name="result")}/phenotype_{phenoname}_{analysis_type}_results.ht'
    if not hfs.exists(f'{all_path}/_SUCCESS'):
        if analysis_type == 'variant':
            fields = ['pop', 'Pvalue']
        else:
            fields = ['pop', 'max_MAF', 'Pvalue', 'Pvalue_Burden', 'Pvalue_SKAT']
        for pop in POPS:
            eur_path = f'{get_aou_saige_results_root(analysis_type=analysis_type, pop="eur", name="result")}/phenotype_{phenoname}/{analysis_type}_results.ht'
            ht = hl.read_table(eur_path)
            ht = ht.annotate(pop = pop)
            ht = ht.select(**{field: ht[field] for field in fields})
            path = f'{get_aou_saige_results_root(analysis_type=analysis_type, pop=pop, name="result")}/phenotype_{phenoname}/{analysis_type}_results.ht'
            if not hfs.exists(f'{path}/_SUCCESS') or pop=='eur':
                continue
            tmp_ht = hl.read_table(path)
            tmp_ht = tmp_ht.annotate(pop = pop)
            tmp_ht = tmp_ht.select(**{field: tmp_ht[field] for field in fields})
            ht = ht.union(tmp_ht, unify=True)
        ht = ht.checkpoint(all_path, _read_if_exists=True)
    ht = hl.read_table(all_path)
    print(ht.count())
    ht.export(all_path.replace('.ht', '.txt.bgz'))
    return all_path.replace('.ht', '.txt.bgz')

def get_phecode_min_p_ht(
    analysis_type: str = "gene",
    result_type: str = "gene",
    test_type: str = "skato",
    pop:str = 'eas',
    filters: bool = True,
    tranche: str = TRANCHE,
):
    """
    Generate min pvalue for each phecode group from variant-level or gene-level test result
    :param str result_type: Type of association test result, `gene` or `variant`
    :param str test_type: Type of gene-level association test result, `skato` or `burden`
    :param bool filters: Whether to use QCed sets
    :param str tranche: Tranche of data to use
    :return: Hail Table of min p-values for each phecode group
    :rtype: Table
    """
    pvalue = P_VALUE_FIELDS[test_type.lower()]
    if filters:
        mt = hl.read_matrix_table(
            get_aou_final_results_path(analysis_type=analysis_type, pop=pop, result_type=result_type, extension='mt')) # TODO: edit this, what filters to add?
    else:
        mt = hl.read_matrix_table(get_aou_final_results_path(analysis_type=analysis_type, pop=pop, result_type=result_type, extension='mt'))

    phecode_info_ht = hl.read_table(PHECODE_PATH)
    mt = mt.annotate_cols(**phecode_info_ht[mt.col_key])

    phecodes = mt.filter_cols(hl.is_defined(mt.category))
    CATEGORIES = phecodes.aggregate_cols(hl.agg.collect_as_set(phecodes.category))
    if result_type == "gene":
        phecodes = phecodes.annotate_rows(
            min_p=hl.agg.group_by(phecodes.category, hl.agg.min(phecodes[pvalue])),
            CHR = phecodes.interval.start.contig,
            POS = phecodes.interval.start.position
        )
        phecodes = phecodes.select_rows(phecodes.CHR, phecodes.POS,
            **{
                f"{category}_min_p": phecodes.min_p.get(category)
                for category in CATEGORIES
            },
        )
    else:
        phecodes = phecodes.annotate_rows(
            CHR=phecodes.locus.contig,
            POS=phecodes.locus.position,
            min_p=hl.agg.group_by(
                phecodes.category, hl.agg.filter(phecodes.SE != 0, hl.agg.min(phecodes.Pvalue))
            )
        )
        phecodes = phecodes.select_rows(phecodes.CHR, phecodes.POS,
            **{
                f"{category}_min_p": phecodes.min_p.get(category)
                for category in CATEGORIES
            }
        )
    return phecodes.rows()

def get_lambda_gc(path, p_field, tag):
    tag = '_Cauchy' if tag == '_cauchy' else tag
    ht = hl.read_table(path)
    if tag != '':
        lambda_gc = hl.eval(ht[f'lambda_gc_maxmaf{tag}'][f'lambda_gc_{p_field.replace("_log10", "")}'])
    else:
        lambda_gc = hl.eval(ht[f'lambda_gc_{p_field.replace("_log10", "")}'])
    return lambda_gc


def main(args):
    ancestries = args.ancestries.split(",")
    analysis_type = "variant" if args.single_variant_only else "gene"
    test_type = 'saige' if analysis_type == 'variant' else 'saige_gene'
    result_type = 'variant' if not args.plot_gene_result else 'gene'
    print(f'Analysis type: {analysis_type}')
    print(f'Result type: {result_type}')
    job_name = 'genome_variant' if analysis_type == 'variant' and result_type !='gene' else 'exome_variant'
    if result_type=='gene' and analysis_type == 'gene':
        job_name = 'gene'
    all_phenos_by_group = read_pickle_dict(f'{DATA_PATH}/phenotype/summary/pheno_by_ancestry_by_group_dict.dict')
    phenos_to_run_by_pop = read_pickle_dict(f'{DATA_PATH}/phenotype/summary/pheno_by_ancestry_dict_{args.sex}.dict')
    phenos_to_run_by_pop_by_group = read_pickle_dict(
        f'{DATA_PATH}/phenotype/summary/pheno_by_ancestry_by_group_dict_{args.sex}.dict')
    phenos_to_run_by_pop_by_group['meta'] = all_phenos_by_group

    backend = hb.ServiceBackend(
        billing_project="all-by-aou", remote_tmpdir=TMP_BUCKET
    )

    b = hb.Batch(
        name=f"aou_{job_name}{'_manhattan' if not args.skip_manhattan_plot else''}_{args.ancestries}",
        backend=backend,
        default_image=args.r_docker_image,
        default_storage="500Mi",
        default_cpu=8,
    )

    if not args.skip_manhattan_plot:
        no_qqplot = 'TRUE' if not args.include_qq_plot else 'FALSE'
        label = 'gene_symbol' if analysis_type=='gene' and result_type == 'gene' else 'locus'
        test_type = 'saige' if analysis_type == 'variant' else 'saige_gene'
        p_method = 'exact' if analysis_type == 'gene' else 'approx_cdf'
        alpha = 6.7e-7 if result_type == 'gene' else 5e-8
        tags = args.tag.split(",") if analysis_type == 'gene' and result_type == 'gene' else ['']
        

        for ancestry in ancestries:
            sex_tag = '' if args.sex == 'both' else f"_{args.sex}"
            software_type = 'gate_' if args.category == 'gate' else ''
            root = f'{ANALYSIS_BUCKET}/{software_type}ht_results{sex_tag}/{ancestry.upper()}'
            category = ''
            if args.category is not None:
                category = args.category
                if category == 'quantitative' or category == 'quantitative_traits':
                    phenos_to_run = set(
                        list(phenos_to_run_by_pop_by_group[ancestry]["lab_measurement"]) +
                        list(phenos_to_run_by_pop_by_group[ancestry]["physical_measurement"]))
                    print(list(phenos_to_run)[0:5])
                elif category == 'gate':
                    phenos_to_run = []
                    for sex_group in ['both', 'male', 'female']:
                        for category in ['onset', 'progression']:
                            phenos_to_run_by_ancestry_by_group = read_pickle_dict(f'{DATA_PATH}/phenotype/summary/pheno_by_ancestry_by_group_dict_{sex_group}.dict')
                            phenos_to_run = phenos_to_run + phenos_to_run_by_ancestry_by_group[ancestry][category]
                    print(list(phenos_to_run)[0:5])
                else: 
                    phenos_to_run = phenos_to_run_by_pop_by_group[ancestry][category]
                    print(phenos_to_run[0:5])
            elif args.phenos is not None:
                phenos_to_run = args.phenos.split(',')
            elif args.pilot:
                phenos_to_run = PILOT_PHENOTYPES
            else:
                phenos_to_run = phenos_to_run_by_pop[ancestry]
                print(f'N phenos: {len(phenos_to_run)}')
                print(phenos_to_run[0:5])

            for i in tqdm(range(len(phenos_to_run))):
                phenoname = list(phenos_to_run)[i]
                for tag in tags:
                    if tag == '_cauchy' and 'SKAT' not in args.p_field: continue
                    input_path = f"{root}/phenotype_{phenoname}/{job_name}_results_{p_method}_expected_p{'_Cauchy' if tag == '_cauchy' else tag}_QC.txt.bgz"
                    print(input_path)

                    if hl.hadoop_exists(input_path):
                        # lambda_gc = hl.eval(hl.read_table(input_path.replace('.txt.bgz', '.ht')).lambda_gc)
                        # print(f'Lambda GC: {lambda_gc}')
                        output_dir = f'{ANALYSIS_BUCKET}/{software_type}{analysis_type}_results/plot/{ancestry.upper()}/phenotype_{phenoname}'
                        print(f'INPUT PATH: {input_path}')
                        print(f'OUTPUT DIR: {output_dir}')

                        if not args.dataproc:
                            input_data = b.read_input(input_path)
                            if args.include_qq_plot:
                                j_lambda = b.new_python_job(name=f'obtain_lambda_gc_{ancestry}_{phenoname}_{analysis_type}')
                                j_lambda.image("hailgenetics/hail:0.2.133-py3.11")
                                j_lambda.memory('lowmem')
                                lambda_gc = j_lambda.call(get_lambda_gc,
                                                        path=input_path.replace(f'_QC.txt.bgz', '.ht').replace('_exact_expected_p_0.001', ''),
                                                        p_field = args.p_field,
                                                        tag = tag
                                    )
                            else:
                                j_lambda = None
                            # b.write_output(lambda_gc.as_str(), f'gs://aou_tmp/lambda_gc/{pop}_{phenoname}_{analysis_type}_{result_type}.txt')

                            j_plot = b.new_job(name=f'create_manhattan_plot_{ancestry}_{phenoname}_{analysis_type}_results')
                            if j_lambda is not None:
                                j_plot.depends_on(j_lambda)
                            j_plot.image(args.r_docker_image)
                            j_plot.command('ls -lha /R')
                            j_plot.memory('highmem')
                            output_files = {}
                            if analysis_type == 'gene' and result_type == 'gene':
                                annotations = ['_pLoF', '_missenseLC', '_synonymous', '']
                                if tag == '_cauchy': annotations = ['']
                                if args.include_qq_plot:
                                    for annt in annotations:
                                        for p_type in ['qqplot', 'qq_and_manhattan']:
                                            output_files[f'{p_type}{annt}.png'] = f'{{root}}{annt}_{args.p_field}_{p_type}.png'
                                else:
                                    for annt in annotations:
                                            output_files[f'manhattan{annt}.png'] = f'{{root}}{annt}_{args.p_field}_manhattan.png'
                            else:
                                if args.include_qq_plot:
                                    for p_type in ['qqplot', 'qq_and_manhattan']:
                                        output_files[f'{p_type}.png'] = f'{{root}}_{args.p_field}_{p_type}.png'
                                else:
                                    output_files={'manhattan.png': f'{{root}}_{args.p_field}_manhattan.png',}
                            j_plot.declare_resource_group(result=output_files)
                            if '_log10' in args.p_field:
                                R_command = f"Rscript /R/manhattan_and_qq_plot.R -a {alpha} -e {args.p_field.replace('_log10', '_expected_log10')} -f {input_data} -g $(<{lambda_gc.as_str()}) -q {no_qqplot} -n {label} -p {args.p_field} -w AC_bin -c CHR -b POS -o {j_plot.result}; "
                                # R_command = f"Rscript /R/manhattan_and_qq_plot.R -a {alpha} -e {args.p_field.replace('_log10', '_expected_log10')} -f {input_data} -g 1 -q {no_qqplot} -n {label} -p {args.p_field} -w AC_bin -c CHR -b POS -o {j_plot.result}; "
                            else:
                                R_command = f"Rscript /R/manhattan_and_qq_plot.R -a {alpha} -f {input_data} -g $(<{lambda_gc.as_str()}) -q {no_qqplot} -n {label} -p {args.p_field} -w AC_bin -c CHR -b POS -o {j_plot.result}; "
                                # R_command = f"Rscript /R/manhattan_and_qq_plot.R -a {alpha} -f {input_data} -g 1 -q {no_qqplot} -n {label} -p {args.p_field} -w AC_bin -c CHR -b POS -o {j_plot.result}; "
                            j_plot.command(R_command)
                            # j_plot.command(f'find io/batch/ -type f')
                            b.write_output(j_plot.result, f'{output_dir}.{ancestry.upper()}.{args.p_field}{tag}.{result_type}')
                        if args.test:
                            break
                    if args.test:
                        break
                if args.test:
                    break

        b.run()

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--single-variant-only", help="Make plot for single variant results", action="store_true"
    )
    parser.add_argument(
        "--plot-gene-result", help="Make plot for gene results", action="store_true"
    )
    parser.add_argument(
        "--meta-analysis", help="Make plot for meta-analysis", action="store_true"
    )
    parser.add_argument(
        "--overwrite", help="Overwrite the files for plotting", action="store_true"
    )
    parser.add_argument(
        "--dataproc", help="Run subsetting on dataproc", action="store_true"
    )
    parser.add_argument(
        "--category",
        help="phenotype category to export",
    )
    parser.add_argument(
        "--result-type", help="result type", default="gene", choices=['gene', 'variant']
    )
    parser.add_argument(
        "--test-type", help="gene test type", default="skato", choices=['skato', 'skat', 'burden']
    )
    parser.add_argument(
        "--ancestries", help="comma-separated list", default="afr,amr,eas,eur,mid,sas"
    )
    parser.add_argument(
        "--tag", help="gene-level results frequency group", default="_0.001,_0.01,_0.0001,_cauchy"
    )
    parser.add_argument(
        "--p-field", help="p value field", default="Pvalue_Burden_log10"
    )
    parser.add_argument(
        "--include-qq-plot", help="include the qqplot component in manhattan plot", action="store_true"
    )
    parser.add_argument(
        "--skip-pop-qq-plot", help="skip running population combined qq plot", action="store_true"
    )
    parser.add_argument(
        "--skip-manhattan-plot", help="skip running manhattan plot", action="store_true"
    )
    parser.add_argument(
        "--skip-phecode-mega-plot", help="skip running phecode mega manhattan plot", action="store_true"
    )
    parser.add_argument(
        "-r-docker-image", help="R docker image to use", action="store_true",
        # default='us-central1-docker.pkg.dev/aou-neale-gwas/plot/aou_man_qq:1.3'
        # default = 'us-central1-docker.pkg.dev/aou-neale-gwas/plot/aou_man_qq:1.6'
        default = 'us-central1-docker.pkg.dev/aou-neale-gwas/plot/aou_man_qq_with_python3.11:latest' # classic manhattan plot
        # default = 'us-central1-docker.pkg.dev/aou-neale-gwas/plot/aou_man_qq_with_python3.11:2.3' # color by AC bin
        # default = 'us-central1-docker.pkg.dev/aou-neale-gwas/plot/aou_man_qq_with_python3.11:2.4' # color by AC bin and stratify into panels
    )
    parser.add_argument(
        "--test", help="Only run a few tests", action="store_true"
    )
    parser.add_argument(
        "--phenos",
        help="Comma-separated list of phenotype names"
    )
    parser.add_argument(
        "--pilot", help="run pilot phenos", action="store_true"
    )
    parser.add_argument(
        "--sex", help="run pilot phenos", default="both"
    )
    args = parser.parse_args()

    main(args)









