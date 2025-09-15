import hail as hl
import hailtop.batch as hb

AC_CUTOFFS = list(range(0, 6)) + [10, 20, 50, 100]
AF_CUTOFFS = sorted([0] + [y * 10 ** x for y in (1, 2, 5) for x in range(-4, 0)] + [0.99])
SIG_THRESHOLD = 5e-8


def format_pheno_dir(pheno):
    return pheno.replace("/", "_")


def get_top_p_from_mt(mt, p, return_ht = True):
    top_p_hit = hl.agg.filter(hl.is_defined(p) & ~hl.is_nan(p),
                              hl.agg.take(mt.entry.annotate(**mt.col), 1, ordering=p))
    mt = mt.annotate_rows(top_p=hl.or_missing(hl.len(top_p_hit) > 0, top_p_hit[0]))
    if return_ht:
        ht = mt.rows()
        return ht.transmute(**ht.top_p)
    else:
        return mt


def get_vep_formatted_data(vep_path: str, legacy_annotations: bool = False):
    from aou_gwas.utils.annotations import annotation_case_builder, annotation_case_builder_ukb_legacy
    from gnomad.utils.vep import process_consequences
    ht = hl.read_table(vep_path)
    ht = process_consequences(ht)
    ht = ht.explode(ht.vep.worst_csq_by_gene_canonical)
    annotation_func = annotation_case_builder_ukb_legacy if legacy_annotations else annotation_case_builder
    return ht.select(
        gene=ht.vep.worst_csq_by_gene_canonical.gene_symbol,
        annotation=annotation_func(ht.vep.worst_csq_by_gene_canonical))



def get_cases_and_controls_from_log(log_format):
    """
    'gs://path/to/result_chr{chrom}_000000001.variant.log'
    """
    cases = controls = -1
    for chrom in range(10, 23):
        try:
            with hfs.open(log_format.format(chrom=chrom)) as f:
                for line in f:
                    line = line.strip()
                    if line.startswith('Analyzing'):
                        fields = line.split()
                        if len(fields) == 6:
                            try:
                                cases = int(fields[1])
                                controls = int(fields[4])
                                break
                            except ValueError:
                                logger.warn(f'Could not load number of cases or controls from {line}.')
                    elif line.endswith('samples were used in fitting the NULL glmm model and are found in sample file') or \
                            line.endswith('samples have been used to fit the glmm null model'):
                        # This is ahead of the case/control count line ("Analyzing ...") above so this should be ok
                        fields = line.split()
                        try:
                            cases = int(fields[0])
                        except ValueError:
                            logger.warn(f'Could not load number of cases or controls from {line}.')
            return cases, controls
        except:
            pass
    return cases, controls


def get_heritability_from_log(log_file, quantitative_trait: bool = False):
    import math
    heritability = -1
    with hfs.open(log_file) as f:
        for line in f:
            if line.startswith('Final'):
                fields = line.strip().split()
                if len(fields) == 4:
                    try:
                        tau = float(fields[2])
                        if quantitative_trait:
                            tau1 = float(fields[1])
                            heritability = tau / (tau1 + tau)
                        else:
                            heritability = tau / (tau + math.pi ** 2 / 3)
                        break
                    except:
                        logger.warn(f'Could not load heritability from {line}.')
    return heritability


def get_saige_version_from_log(null_glmm_log):
    version = 'NA'
    with hfs.open(null_glmm_log) as f:
        for line in f:
            if line.startswith('other attached packages:'):
                try:
                    line2 = f.readline()
                    packages = line2.strip().split()
                    version = [x for x in packages if 'SAIGE' in x][0]
                except:
                    logger.warning(f'Could not load version number from {line2} in {null_glmm_log}.')
    return version


def get_inverse_normalize_status(null_glmm_log):
    status = 'Unknown'
    with hfs.open(null_glmm_log) as f:
        for line in f:
            if line.startswith('$invNormalize'):
                try:
                    status = f.readline().strip().split()[1]
                except:
                    logger.warning(f'Could not load inv_norm status from {line} in {null_glmm_log}.')
    return status.capitalize()


def get_saige_timing_grep(all_files):
    try:
        grep_results = hl.grep('Analysis took', all_files, max_count=int(1e8), show=False)
    except hl.utils.java.FatalError:
        return
    if sum([len(x) for x in grep_results.values()]) > 5e7:
        logger.warning(f'Got more than 5e7 values in {all_files[0]}, etc. Check this!')
    for log, result in grep_results.items():
        try:
            timing = float(result[0].split()[2])
        except:
            logger.warning(f'Could not load timing from {result} in {log}.')
            continue
        chrom, pos = log.rsplit('.', 2)[0].rsplit('_', 2)[1:3]
        yield f'{chrom}:{pos}', timing


def get_null_model_timing(null_glmm_log):
    cpu = wall = 'NA'
    with hfs.open(null_glmm_log) as f:
        for line in f:
            if line.startswith('t_end - t_begin'):
                try:
                    f.readline()
                    line2 = f.readline()
                    cpu, _, wall = line2.strip().split()
                except:
                    logger.warning(f'Could not load null model timings from {line2} in {null_glmm_log}.')
    return cpu, wall


def union_mts_by_tree(all_mts, temp_dir, debug=False):
    chunk_size = int(len(all_mts) ** 0.5) + 1
    outer_mts = []
    for i in range(chunk_size):
        if i * chunk_size >= len(all_mts): break
        mt = all_mts[i * chunk_size]
        for j in range(1, chunk_size):
            if i * chunk_size + j >= len(all_mts): break
            try:
                mt = mt.union_cols(all_mts[i * chunk_size + j], row_join_type='outer')
            except:
                if debug:
                    print(f'problem with {i * chunk_size} and {i * chunk_size + j}')
                    mt.describe()
                    all_mts[i * chunk_size + j].describe()
                raise
        outer_mts.append(mt.checkpoint(f'{temp_dir}/temp_output_{i}.mt', overwrite=True))
    mt = outer_mts[0]
    for next_mt in outer_mts[1:]:
        mt = mt.union_cols(next_mt, row_join_type='outer')
    return mt


def union_hts_by_tree(all_hts, temp_dir, debug=False, inner_mode = 'overwrite'):
    chunk_size = int(len(all_hts) ** 0.5) + 1
    outer_hts = []
    for i in range(chunk_size):
        if i * chunk_size >= len(all_hts): break
        hts = all_hts[i * chunk_size:(i + 1) * chunk_size]
        try:
            if isinstance(hts[0], str):
                hts = list(map(lambda x: hl.read_table(x), hts))
            ht = hts[0].union(*hts[1:], unify=True)
        except:
            if debug:
                print(f'problem in range {i * chunk_size}-{i * chunk_size + chunk_size}')
                _ = [ht.describe() for ht in hts]
            raise
        outer_hts.append(ht.checkpoint(f'{temp_dir}/temp_output_{i}.ht', **{inner_mode: True}))
    return outer_hts[0].union(*outer_hts[1:], unify=True)


def get_files_in_parent_directory(parent_dir, fname: str = 'variant_results.ht'):
    all_outputs = []
    for directory in parent_dir:
        if not directory['is_dir']:
            continue
        file_path = f'{directory["path"]}/{fname}'
        if hfs.exists(f'{file_path}/_SUCCESS'):
            all_outputs.append(file_path)
    return all_outputs


def union_ht(all_hts, col_fields, pheno_dict, temp_dir, inner_mode: str = 'overwrite'):
    print(f'Unioning {len(all_hts)} HTs...')
    ht = union_hts_by_tree(all_hts, temp_dir, inner_mode=inner_mode)
    return ht.annotate(**pheno_dict[ht.key.select(*col_fields)])


def pull_out_col_keys(all_hts, row_keys, col_keys):
    rekeyed_hts = []
    for ht in all_hts:
        ht2 = ht.head(1)
        glob = ht2.aggregate(hl.agg.take(hl.struct(**{x: ht2[x] for x in col_keys}), 1)[0], _localize=False)
        rekeyed_hts.append(ht.key_by(*row_keys).drop(*col_keys).annotate_globals(**glob))
    return rekeyed_hts


def join_pheno_hts_to_mt(all_hts, row_keys, col_keys, temp_dir = None, inner_mode: str = 'overwrite',
                         repartition_final: int = None):
    rekeyed_hts = pull_out_col_keys(all_hts, row_keys, col_keys)
    # mt = mwzj_hts_by_tree_on_batch(batch, rekeyed_hts, temp_dir, col_keys, debug=True,
    #                               inner_mode=inner_mode, repartition_final=repartition_final)
    mt = mwzj_hts_by_tree(rekeyed_hts, temp_dir, col_keys, debug=True,
                                  inner_mode=inner_mode, repartition_final=repartition_final)
    print(f'Unioned MTs...')
    return mt


def unify_saige_ht_schema(ht, result_type, patch_case_control_count: str = ''):
    """

    :param Table ht:
    :param str patch_case_control_count: Path to file (hack to get cases and controls back if loading later)
    :return:
    :rtype: Table
    """
    # assert ht.head(1).annotation.collect()[0] is None, f'failed at {patch_case_control_count}'
    if result_type == 'variant':
        ht = ht.select(**{'MarkerID': ht['MarkerID'], 'AC_Allele2': ht['AC_Allele2'], 'AF_Allele2': ht['AF_Allele2'], 'MissingRate': ht['MissingRate'], 'BETA':hl.float64(ht.BETA),
            'SE': hl.float64(ht.SE), 'var': ht.var,'p.value.NA':ht['p.value.NA'], 'Is.SPA': ht['Is.SPA'], 'AF_case': ht['AF_case'], 'AF_ctrl': ht['AF_ctrl'], 'Pvalue': ht['Pvalue'],
        })
    else:
        ht.describe()
        ht = ht.select('Pvalue', 'Pvalue_Burden', 'Pvalue_SKAT', 'Pvalue_log10', 'Pvalue_Burden_log10', 'Pvalue_SKAT_log10', 'CHR', 'POS',
              'BETA_Burden', 'SE_Burden','MAC',  'Number_rare', 'Number_ultra_rare', 'total_variants', 'interval')


    # ht2 = ht.head(1)
    # pheno_key_dict = dict(ht2.aggregate(hl.agg.take(ht2.key, 1)[0]))
    # if patch_case_control_count:
    #     if not ht.n_cases.collect()[0]:
    #         directory, tpc, _ = patch_case_control_count.rsplit('/', 2)
    #
    #         pheno_results_dir = get_pheno_output_path(directory, pheno_key_dict, '', legacy=True)
    #         prefix = get_results_prefix(pheno_results_dir, pheno_key_dict, '{chrom}', 1, legacy=True)
    #         saige_log = f'{prefix}.variant.log'
    #         cases, controls = get_cases_and_controls_from_log(saige_log)
    #         print(f'Patched pheno: {tpc}. Got {cases} cases and {controls} controls.')
    #         if cases == -1: cases = hl.null(hl.tint)
    #         if controls == -1: controls = hl.null(hl.tint)
    #         ht = ht.annotate_globals(n_cases=cases, n_controls=controls)
    if 'heritability' not in list(ht.globals):
        ht = ht.annotate_globals(heritability=hl.null(hl.tfloat64))
    if 'saige_version' not in list(ht.globals):
        ht = ht.annotate_globals(saige_version=hl.null(hl.tstr))
    if 'log_pvalue' not in list(ht.globals):
        ht = ht.annotate_globals(log_pvalue=False)
    return ht


def stringify_pheno_key_dict(pheno_key_dict, format_phenocode_field: bool = False, delimiter='-'):
    return delimiter.join([format_pheno_dir(pheno_key_dict[x])
                           if x == 'phenocode' and format_phenocode_field
                           else pheno_key_dict[x] for x in PHENO_KEY_FIELDS if x in pheno_key_dict])



def get_results_prefix(pheno_results_dir, pheno_key_dict, chromosome, start_pos, legacy: bool = False):
    prefix = f'{pheno_results_dir}/result_'
    # prefix = f'{pheno_results_dir}/result_megabase_span_'  # TODO: add for megabase spanning bug
    if legacy:
        prefix += format_pheno_dir(pheno_key_dict["phenocode"])
    else:
        prefix += stringify_pheno_key_dict(pheno_key_dict, True)
    return f'{prefix}_{chromosome}_{str(start_pos).zfill(9)}'


def get_pheno_output_path(pheno_export_dir, pheno_coding_trait, extension = '.tsv', legacy: bool = False):
    if legacy:
        extended_suffix = pheno_coding_trait['coding']
    else:
        extended_suffix = f'{pheno_coding_trait["pheno_sex"]}-{pheno_coding_trait["coding"]}-{pheno_coding_trait["modifier"]}'
    return f'{pheno_export_dir}/{pheno_coding_trait["trait_type"]}-{format_pheno_dir(pheno_coding_trait["phenocode"])}-{extended_suffix}{extension}'


def recode_pkd_to_legacy(pheno_key_dict_list):
    for pheno_key_dict in pheno_key_dict_list:
        recode_single_pkd_to_legacy(pheno_key_dict)
    return pheno_key_dict_list


def recode_single_pkd_to_legacy(pheno_key_dict):
    if pheno_key_dict['trait_type'] == 'icd10':
        pheno_key_dict['trait_type'] = 'icd_all'
        pheno_key_dict['coding'] = 'icd10'
    elif pheno_key_dict['trait_type'] == 'phecode':
        pheno_key_dict['coding'] = pheno_key_dict['pheno_sex']
    elif pheno_key_dict['trait_type'] == 'biomarkers':
        pheno_key_dict['coding'] = pheno_key_dict['phenocode']
    else:
        if pheno_key_dict['phenocode'] == 'whr':
            pheno_key_dict['coding'] = 'whr'
        else:
            pheno_key_dict['coding'] = pheno_key_dict['coding'] if pheno_key_dict['coding'] else pheno_key_dict['modifier']
    del pheno_key_dict['pheno_sex']
    del pheno_key_dict['modifier']


def recode_pkd_to_new(pheno_key_dict_list):
    for pheno_key_dict in pheno_key_dict_list:
        recode_single_pkd_to_new(pheno_key_dict)
    return pheno_key_dict_list


def recode_single_pkd_to_new(pheno_key_dict):
    new_dict = {}
    SEX = 'both_sexes'
    if pheno_key_dict['trait_type'] == 'icd_all':
        new_dict['trait_type'] = 'icd10'
        new_dict['phenocode'] = pheno_key_dict['phenocode']
        new_dict['pheno_sex'] = SEX
        new_dict['coding'] = ''
        new_dict['modifier'] = ''
    else:
        new_dict['trait_type'] = pheno_key_dict['trait_type']
        new_dict['phenocode'] = pheno_key_dict['phenocode']
        if pheno_key_dict['trait_type'] == 'phecode':
            new_dict['pheno_sex'] = pheno_key_dict['coding']
            new_dict['coding'] = ''
            new_dict['modifier'] = ''
        else:
            new_dict['pheno_sex'] = SEX
            if pheno_key_dict['trait_type'] == 'categorical':
                new_dict['coding'] = pheno_key_dict['coding']
                new_dict['modifier'] = ''
            elif pheno_key_dict['trait_type'] == 'continuous':
                if pheno_key_dict['phenocode'] == 'whr':
                    new_dict['coding'] = ''
                    new_dict['modifier'] = 'irnt'
                else:
                    new_dict['coding'] = ''
                    new_dict['modifier'] = pheno_key_dict['coding']
            else:
                new_dict['coding'] = ''
                new_dict['modifier'] = ''
    return new_dict


def unify_saige_ht_variant_schema(ht):
    shared = ("MarkerID", "AC_Allele2", "AF_Allele2", "MissingRate", "BETA", "SE", "var", "p.value.NA", "Is.SPA", "AF_case", "AF_ctrl", "Pvalue")
    ht = ht.select(*shared)
    if 'CHR' in list(ht.row):
        ht = ht.drop("CHR", "POS", "Allele1", "Allele2")
    # if 'N' in list(ht.row):
    #     shared = ('MarkerID', 'AC', 'AF', 'N', 'BETA', 'SE', 'Tstat', 'varT', 'varTstar')
    # else:
    #     shared = ('MarkerID', 'AC', 'AF', 'BETA', 'SE')
    # new_floats = ('AF_case', 'AF_ctrl')
    # new_ints = ('N.Cases', 'N.Controls')
    # shared_end = ('Pvalue')
    # if 'AF.Cases' not in list(ht.row):
    #     # TODO: store AF in AF.cases
    #     ht = ht.select(*shared, **{field: hl.null(hl.tfloat64) for field in new_floats},
    #                    **{field: hl.null(hl.tint32) for field in new_ints},
    #                    **{field: ht[field] for field in shared_end})
    # else:
    #     ht = ht.select(*shared, *new_floats, *new_ints, *shared_end)
    return ht.annotate(SE=hl.float64(ht.SE))


def unify_saige_burden_ht_schema(ht, pop='other'):
    ht.describe()
    if 'phenoname' not in list(ht.row_value):
        phenoname = hl.eval(ht.phenoname)
        ht = ht.drop('phenoname')
        ht = ht.annotate(phenoname = phenoname)
        ht = ht.key_by('gene_id', 'gene_symbol', 'annotation', 'max_MAF', 'phenoname')
    if pop != 'meta':
        shared = ( 'Pvalue', 'Pvalue_Burden', 'Pvalue_SKAT', 'Pvalue_log10', 'Pvalue_Burden_log10', 'Pvalue_SKAT_log10', 'CHR', 'POS',
                  'BETA_Burden', 'SE_Burden','MAC',  'Number_rare', 'Number_ultra_rare', 'total_variants', 'interval')
        new_ints = ('MAC_case', 'MAC_control')
        if 'MAC_case' not in list(ht.row):
            ht = ht.select(*shared, **{field: hl.missing(hl.tint32) for field in new_ints})
        else:
            ht = ht.select(*shared, *new_ints)
    else:
        shared = ('Pvalue', 'Pvalue_Burden', 'Pvalue_SKAT', 'Pvalue_log10', 'Pvalue_Burden_log10', 'Pvalue_SKAT_log10',
                  'META_Stats_SKATO', 'META_Stats_SKAT', 'META_Stats_Burden')
        new_ints = ('CHR', 'POS')
        if 'interval' not in list(ht.row_value):
            ht = ht.select(*shared, CHR=hl.missing(hl.tstr), POS = hl.missing(hl.tint32))
        else:
            ht = ht.select(*shared, *new_ints)
    return ht


def get_n_even_intervals(n):
    ref = hl.default_reference()
    genome_size = sum(ref.lengths.values())
    partition_size = int(genome_size / n) + 1
    global_locus_intervals = [
        (x * partition_size, min(x * partition_size + partition_size, genome_size - 1))
        for x in range(n)
    ]
    return list((
        map(
            lambda x: hl.Interval(
                hl.get_reference('default').locus_from_global_position(x[0]),
                hl.get_reference('default').locus_from_global_position(x[1])
            ), global_locus_intervals
        )
    ))




def mwzj_hts_by_tree(all_hts, temp_dir, globals_for_col_key, debug=False, inner_mode = 'overwrite',
                     repartition_final: int = None):
    chunk_size = int(len(all_hts) ** 0.5) + 1
    outer_hts = []

    checkpoint_kwargs = {inner_mode: True}
    if repartition_final is not None:
        intervals = get_n_even_intervals(repartition_final)
        checkpoint_kwargs['_intervals'] = intervals

    if debug: print(f'Running chunk size {chunk_size}...')
    for i in range(chunk_size):
        if i * chunk_size >= len(all_hts): break
        hts = all_hts[i * chunk_size:(i + 1) * chunk_size]
        if debug: print(f'Going from {i * chunk_size} to {(i + 1) * chunk_size} ({len(hts)} HTs)...')
        try:
            if isinstance(hts[0], str):
                hts = list(map(lambda x: hl.read_table(x), hts))
            ht = hl.Table.multi_way_zip_join(hts, 'row_field_name', 'global_field_name')
        except:
            if debug:
                print(f'problem in range {i * chunk_size}-{i * chunk_size + chunk_size}')
                _ = [ht.describe() for ht in hts]
            raise
        outer_hts.append(ht.checkpoint(f'{temp_dir}/temp_output_{i}.ht', **checkpoint_kwargs))
    ht = hl.Table.multi_way_zip_join(outer_hts, 'row_field_name_outer', 'global_field_name_outer')
    ht = ht.transmute(inner_row=hl.flatmap(lambda i: hl.coalesce(ht.row_field_name_outer[i].row_field_name, hl.range(0, hl.len(ht.global_field_name_outer[i].global_field_name))
                                                   .map(lambda _: hl.null(ht.row_field_name_outer[i].row_field_name.dtype.element_type))),
                                           # hl.cond(hl.is_missing(ht.row_field_name_outer[i].row_field_name),
                                           #         hl.range(0, hl.len(ht.global_field_name_outer[i].global_field_name))
                                           #         .map(lambda _: hl.null(ht.row_field_name_outer[i].row_field_name.dtype.element_type)),
                                           #         ht.row_field_name_outer[i].row_field_name),
                                           hl.range(hl.len(ht.global_field_name_outer))))
    ht = ht.transmute_globals(inner_global=hl.flatmap(lambda x: x.global_field_name, ht.global_field_name_outer))
    mt = ht._unlocalize_entries('inner_row', 'inner_global', globals_for_col_key)
    return mt


def generate_lambda_ht_by_freq(mt):
    af_cases = mt['AF.Cases']
    ac_cases = af_cases * mt.n_cases * 2
    af_total = mt['AF_Allele2']
    p_value_field = hl.exp(mt.Pvalue)
    breakdown_dict_tuple = (('by_case', {'ac': ac_cases}), ('by', {'af': af_total}))
    mt = mt.annotate_cols(sumstats_qc=generate_qc_lambda_aggregator(breakdown_dict_tuple, p_value_field)).annotate_globals(ac_cutoffs=AC_CUTOFFS, af_cutoffs=AF_CUTOFFS)
    return mt.cols()


def generate_qc_lambda_aggregator(breakdown_dict_tuple, p_value_field):
    return hl.struct(**{
        f'{metric}_{breakdown}_{flavor}': [hl.agg.filter(breakdown_dict[flavor] >= cutoff, agg) for cutoff in cutoffs]
        for flavor, cutoffs in (('ac', AC_CUTOFFS), ('af', AF_CUTOFFS))
        for breakdown, breakdown_dict in breakdown_dict_tuple if flavor in breakdown_dict
        for metric, agg in (
            ('lambda_gc', hl.methods.statgen._lambda_gc_agg(p_value_field)),
            ('n_variants', hl.agg.count()),
            ('n_sig', hl.agg.count_where(p_value_field < SIG_THRESHOLD))
        )
    })


def explode_lambda_ht(ht, by='ac'):
    ac_ht = ht.annotate(sumstats_qc=ht.sumstats_qc.select(*[x for x in ht.sumstats_qc.keys() if f'_{by}' in x]))
    ac_ht = ac_ht.annotate(index_ac=hl.zip_with_index(ac_ht[f'{by}_cutoffs'])).explode('index_ac')
    ac_ht = ac_ht.transmute(**{by: ac_ht.index_ac[1]},
                            **{x: ac_ht.sumstats_qc[x][ac_ht.index_ac[0]] for x in ac_ht.sumstats_qc})
    return ac_ht

def pull_out_fields_from_entries(mt, shared_fields, index='rows', agg_funcs=None):
    func = mt.annotate_rows if index == 'rows' else mt.annotate_cols
    if agg_funcs is None:
        agg_funcs = [hl.agg.take for _ in shared_fields]
    elif not isinstance(agg_funcs, list):
        agg_funcs = [agg_funcs for _ in shared_fields]
    mt = func(**{f'_{field}': agg_func(mt[field]) for agg_func, field in zip(agg_funcs, shared_fields)})
    return mt.drop(*shared_fields).rename({f'_{field}': field for field in shared_fields})

def mwzj_hts_by_tree_on_batch(batch: hb.Batch, all_hts, temp_dir, globals_for_col_key, debug=False, inner_mode = 'overwrite',
                     repartition_final: int = None):
    def _func1(chunk_size, outer_hts, debug):
        hts = all_hts[i * chunk_size:(i + 1) * chunk_size]
        if debug: print(f'Going from {i * chunk_size} to {(i + 1) * chunk_size} ({len(hts)} HTs)...')
        try:
            if isinstance(hts[0], str):
                hts = list(map(lambda x: hl.read_table(x), hts))
            ht = hl.Table.multi_way_zip_join(hts, 'row_field_name', 'global_field_name')
        except:
            if debug:
                print(f'problem in range {i * chunk_size}-{i * chunk_size + chunk_size}')
                _ = [ht.describe() for ht in hts]
            raise
        outer_hts.append(ht.checkpoint(f'{temp_dir}/temp_output_{i}.ht', **checkpoint_kwargs))
        return outer_hts

    def _func2(outer_hts, ):
        ht = hl.Table.multi_way_zip_join(outer_hts, 'row_field_name_outer', 'global_field_name_outer')
        ht = ht.transmute(inner_row=hl.flatmap(lambda i:
                                               hl.cond(hl.is_missing(ht.row_field_name_outer[i].row_field_name),
                                                       hl.range(0,
                                                                hl.len(ht.global_field_name_outer[i].global_field_name))
                                                       .map(lambda _: hl.null(ht.row_field_name_outer[
                                                                                  i].row_field_name.dtype.element_type)),
                                                       ht.row_field_name_outer[i].row_field_name),
                                               hl.range(hl.len(ht.global_field_name_outer))))
        ht = ht.transmute_globals(inner_global=hl.flatmap(lambda x: x.global_field_name, ht.global_field_name_outer))
        mt = ht._unlocalize_entries('inner_row', 'inner_global', globals_for_col_key)
        return mt
    chunk_size = int(len(all_hts) ** 0.5) + 1
    outer_hts = []

    checkpoint_kwargs = {inner_mode: True}
    if repartition_final is not None:
        intervals = get_n_even_intervals(repartition_final)
        checkpoint_kwargs['_intervals'] = intervals

    if debug: print(f'Running chunk size {chunk_size}...')
    job_lst = []
    for i in range(chunk_size):
        if i * chunk_size >= len(all_hts): break
        j = batch.new_python_job(name = f'job_{i}')
        outer_hts = j.call(_func1,
                           chunk_size=chunk_size,
                           outer_hts=outer_hts,
                           debug=debug)
        job_lst.append(j)
    j_final = batch.new_python_job(name='final_merge')
    mt = j_final.call(_func2,
                      outer_hts=outer_hts)
    j_final.depends_on(*job_lst)

    return mt

def all_and_leave_one_out(x, pop_array, all_f=hl.sum, loo_f=lambda i, x: hl.sum(x) - hl.or_else(x[i], 0)):
    """
    Applies a function to an input array for all populations, and for each of leave-one-out populations.

    :param x: Input array
    :param pop_array: Population array
    :param all_f: Function for all populations. It takes the input array and returns a new value
    :param loo_f: Function for each of leave-one-out populations. It takes an index of leave-one-out
                  population and the input array, and returns an array of new values.
    ...
    :return: Array of new values for all populations and for each of leave-one-out populations.
    :rtype: ArrayExpression
    """
    arr = hl.array([all_f(x)])
    arr = arr.extend(hl.map(lambda i: loo_f(i, x), hl.range(hl.len(pop_array))))
    return hl.or_missing(hl.any(hl.is_defined, x), arr)


def run_meta_analysis(mt: hl.MatrixTable,
                      beta_field: str = 'BETA',
                      se_field: str = 'SE',
                      af_allele2_field: str = 'AF_Allele2',
                      af_case_field: str='AF_case',
                      af_ctrl_field: str='AF_ctrl'):
    """
    Run inverse-variance fixed-effect meta-analysis for a given MatrixTable.

    :param mt: Input MatrixTable, formatted similar to `load_final_sumstats_mt()`
    ...
    :return: Result MatrixTable with `meta_analysis` entries and `meta_analysis_data` columns.
    :rtype: MatrixTable
    """
    # Annotate per-entry sample size
    def get_n(pheno_data, i):
        return pheno_data[i].n_cases + hl.or_else(pheno_data[i].n_controls, 0)

    mt = mt.annotate_entries(
        summary_stats=hl.map(
            lambda x: x[1].annotate(N=hl.or_missing(hl.is_defined(x[1]), get_n(mt.pheno_data, x[0]))),
            hl.enumerate(mt.summary_stats),
        )
    )

    # Run fixed-effect meta-analysis (all + leave-one-out)
    mt = mt.annotate_entries(
        unnorm_beta=mt.summary_stats[beta_field] / (mt.summary_stats[se_field] ** 2), inv_se2=1 / (mt.summary_stats[se_field] ** 2)
    )
    mt = mt.annotate_entries(
        sum_unnorm_beta=all_and_leave_one_out(mt.unnorm_beta, mt.pheno_data.pop),
        sum_inv_se2=all_and_leave_one_out(mt.inv_se2, mt.pheno_data.pop),
    )
    mt = mt.transmute_entries(
        META_BETA=mt.sum_unnorm_beta / mt.sum_inv_se2, META_SE=hl.map(lambda x: hl.sqrt(1 / x), mt.sum_inv_se2)
    )
    mt = mt.annotate_entries(
        META_Pvalue=hl.map(lambda x: hl.log(2) + hl.pnorm(x, log_p=True), -hl.abs(mt.META_BETA / mt.META_SE))
    )

    # Run heterogeneity test (Cochran's Q)
    mt = mt.annotate_entries(
        META_Q=hl.map(lambda x: hl.sum((mt.summary_stats[beta_field] - x) ** 2 * mt.inv_se2), mt.META_BETA),
        variant_exists=hl.map(lambda x: ~hl.is_missing(x), mt.summary_stats[beta_field]),
    )
    mt = mt.annotate_entries(META_N_pops=all_and_leave_one_out(mt.variant_exists, mt.pheno_data.pop))
    # filter Q-values when N_pops == 1
    mt = mt.annotate_entries(
        META_Q=hl.map(lambda i: hl.or_missing(mt.META_N_pops[i] > 1, mt.META_Q[i]), hl.range(hl.len(mt.META_Q)))
    )
    mt = mt.annotate_entries(
        META_Pvalue_het=hl.map(
            lambda i: hl.pchisqtail(mt.META_Q[i], mt.META_N_pops[i] - 1, log_p=True), hl.range(hl.len(mt.META_Q))
        )
    )

    # Add other annotations
    if af_case_field is not None or af_ctrl_field is not None:
        mt = mt.annotate_entries(
            ac_cases=hl.map(lambda x: x[af_case_field] * x.N, mt.summary_stats),
            ac_controls=hl.map(lambda x: x[af_ctrl_field] * x.N, mt.summary_stats),
            META_AC_Allele2=all_and_leave_one_out(mt.summary_stats[af_allele2_field] * mt.summary_stats.N, mt.pheno_data.pop),
            META_N=all_and_leave_one_out(mt.summary_stats.N, mt.pheno_data.pop),
        )


    mt = mt.annotate_entries(
        META_AF_Allele2=mt.META_AC_Allele2 / mt.META_N,
        META_AF_Cases=all_and_leave_one_out(mt.ac_cases, mt.pheno_data.pop) / mt.META_N,
        META_AF_Controls=all_and_leave_one_out(mt.ac_controls, mt.pheno_data.pop) / mt.META_N,
    )
    mt = mt.drop(
        "unnorm_beta", "inv_se2", "variant_exists", "ac_cases", "ac_controls", "summary_stats", "META_AC_Allele2"
    )

    # Format everything into array<struct>
    def is_finite_or_missing(x):
        return hl.or_missing(hl.is_finite(x), x)

    meta_fields = ["BETA", "SE", "Pvalue", "Q", "Pvalue_het", "N", "N_pops", "AF_Allele2", "AF_Cases", "AF_Controls"]
    mt = mt.transmute_entries(
        meta_analysis=hl.map(
            lambda i: hl.struct(**{field: is_finite_or_missing(mt[f"META_{field}"][i]) for field in meta_fields}),
            hl.range(hl.len(mt.META_BETA)),
        )
    )

    col_fields = ["n_cases", "n_controls"]
    mt = mt.annotate_cols(
        **{field: all_and_leave_one_out(mt.pheno_data[field], mt.pheno_data.pop) for field in col_fields}
    )
    col_fields += ["pop"]
    mt = mt.annotate_cols(
        pop=all_and_leave_one_out(
            mt.pheno_data.pop,
            mt.pheno_data.pop,
            all_f=lambda x: x,
            loo_f=lambda i, x: hl.filter(lambda y: y != x[i], x),
        )
    )
    mt = mt.transmute_cols(
        meta_analysis_data=hl.map(
            lambda i: hl.struct(**{field: mt[field][i] for field in col_fields}), hl.range(hl.len(mt.pop))
        )
    )

    return mt


def run_gene_meta_analysis(mt: hl.MatrixTable,
                           beta_field: str = 'BETA',
                           se_field: str = 'SE',
                           af_allele2_field: str = 'AF_Allele2',
                           ac_case_field: str='MAC_case',
                           ac_ctrl_field: str='MAC_ctrl'):
    """
    Run inverse-variance fixed-effect meta-analysis for a given MatrixTable.

    :param mt: Input MatrixTable, formatted similar to `load_final_sumstats_mt()`
    ...
    :return: Result MatrixTable with `meta_analysis` entries and `meta_analysis_data` columns.
    :rtype: MatrixTable
    """
    # Annotate per-entry sample size
    def get_n(pheno_data, i):
        return pheno_data[i].n_cases + hl.or_else(pheno_data[i].n_controls, 0)

    mt = mt.annotate_entries(
        summary_stats=hl.map(
            lambda x: x[1].annotate(N=hl.or_missing(hl.is_defined(x[1]), get_n(mt.pheno_data, x[0]))),
            hl.enumerate(mt.summary_stats),
        )
    )

    # Run fixed-effect meta-analysis (all + leave-one-out)
    mt = mt.annotate_entries(
        unnorm_beta=mt.summary_stats[beta_field] / (mt.summary_stats[se_field] ** 2), inv_se2=1 / (mt.summary_stats[se_field] ** 2)
    )
    mt = mt.annotate_entries(
        sum_unnorm_beta=all_and_leave_one_out(mt.unnorm_beta, mt.pheno_data.pop),
        sum_inv_se2=all_and_leave_one_out(mt.inv_se2, mt.pheno_data.pop),
    )
    mt = mt.transmute_entries(
        META_BETA=mt.sum_unnorm_beta / mt.sum_inv_se2, META_SE=hl.map(lambda x: hl.sqrt(1 / x), mt.sum_inv_se2)
    )
    mt.META_BETA.show()
    mt.META_SE.show()
    mt = mt.annotate_entries(
        META_Pvalue=hl.map(lambda x: hl.log(2) + hl.pnorm(x, log_p=True), -hl.abs(mt.META_BETA / mt.META_SE))
    )

    # Run heterogeneity test (Cochran's Q)
    mt = mt.annotate_entries(
        META_Q=hl.map(lambda x: hl.sum((mt.summary_stats[beta_field] - x) ** 2 * mt.inv_se2), mt.META_BETA),
        variant_exists=hl.map(lambda x: ~hl.is_missing(x), mt.summary_stats[beta_field]),
    )
    mt = mt.annotate_entries(META_N_pops=all_and_leave_one_out(mt.variant_exists, mt.pheno_data.pop))
    # filter Q-values when N_pops == 1
    mt = mt.annotate_entries(
        META_Q=hl.map(lambda i: hl.or_missing(mt.META_N_pops[i] > 1, mt.META_Q[i]), hl.range(hl.len(mt.META_Q)))
    )
    mt = mt.annotate_entries(
        META_Pvalue_het=hl.map(
            lambda i: hl.pchisqtail(mt.META_Q[i], mt.META_N_pops[i] - 1, log_p=True), hl.range(hl.len(mt.META_Q))
        ),
        META_N=all_and_leave_one_out(mt.summary_stats.N, mt.pheno_data.pop),
    )

    # Add other annotations
    # if af_case_field is not None or af_ctrl_field is not None:
    #     mt = mt.annotate_entries(
    #         ac_cases=hl.map(lambda x: x[af_case_field] * x.N, mt.summary_stats),
    #         ac_controls=hl.map(lambda x: x[af_ctrl_field] * x.N, mt.summary_stats),
    #         META_AC_Allele2=all_and_leave_one_out(mt.summary_stats[af_allele2_field] * mt.summary_stats.N, mt.pheno_data.pop),
    #         META_N=all_and_leave_one_out(mt.summary_stats.N, mt.pheno_data.pop),
    #     )
    #
    #
    # mt = mt.annotate_entries(
    #     META_AF_Allele2=mt.META_AC_Allele2 / mt.META_N,
    #     META_AF_Cases=all_and_leave_one_out(mt.ac_cases, mt.pheno_data.pop) / mt.META_N,
    #     META_AF_Controls=all_and_leave_one_out(mt.ac_controls, mt.pheno_data.pop) / mt.META_N,
    # )
    # mt = mt.drop(
    #     "unnorm_beta", "inv_se2", "variant_exists", "ac_cases", "ac_controls", "summary_stats", "META_AC_Allele2"
    # )
    mt = mt.drop(
        "unnorm_beta", "inv_se2", "variant_exists",  "summary_stats",
    )
    #
    # Format everything into array<struct>
    def is_finite_or_missing(x):
        return hl.or_missing(hl.is_finite(x), x)

    # meta_fields = ["BETA", "SE", "Pvalue", "Q", "Pvalue_het", "N", "N_pops", "AF_Allele2", "AF_Cases", "AF_Controls"]
    meta_fields = ["BETA", "SE", "Pvalue", "Q", "Pvalue_het", "N", "N_pops"]
    mt = mt.transmute_entries(
        meta_analysis=hl.map(
            lambda i: hl.struct(**{field: is_finite_or_missing(mt[f"META_{field}"][i]) for field in meta_fields}),
            hl.range(hl.len(mt.META_BETA)),
        )
    )

    col_fields = ["n_cases", "n_controls"]
    mt = mt.annotate_cols(
        **{field: all_and_leave_one_out(mt.pheno_data[field], mt.pheno_data.pop) for field in col_fields}
    )
    col_fields += ["pop"]
    mt = mt.annotate_cols(
        pop=all_and_leave_one_out(
            mt.pheno_data.pop,
            mt.pheno_data.pop,
            all_f=lambda x: x,
            loo_f=lambda i, x: hl.filter(lambda y: y != x[i], x),
        )
    )
    mt = mt.transmute_cols(
        meta_analysis_data=hl.map(
            lambda i: hl.struct(**{field: mt[field][i] for field in col_fields}), hl.range(hl.len(mt.pop))
        )
    )

    return mt



def run_stouffer(mt):
    P_FIELDS = ['Pvalue', 'Pvalue_Burden', 'Pvalue_SKAT']
    P_TESTS = {'SKATO': 'Pvalue', 'Burden': 'Pvalue_Burden', 'SKAT':'Pvalue_SKAT'}
    def _edit_pvalue(p):
        return hl.if_else(p> 0.99, 0.99, p)

    # Annotate per-entry sample size
    def get_n(pheno_data, i):
        return pheno_data[i].n_cases + hl.or_else(pheno_data[i].n_controls, 0)

    def get_n_eff(pheno_data, i):
        n_controls = hl.or_else(pheno_data[i].n_controls, 0)
        n_cases = pheno_data[i].n_cases
        return (4 * n_cases * n_controls) / (n_cases + n_controls)



    mt = mt.annotate_entries(
        summary_stats=hl.map(
            lambda x: x[1].annotate(N=hl.or_missing(hl.is_defined(x[1]), get_n(mt.pheno_data, x[0]))),
            hl.enumerate(mt.summary_stats),
        )
    )

    mt = mt.annotate_entries(
        summary_stats=hl.map(
            lambda x: x[1].annotate(N_eff=hl.or_missing(hl.is_defined(x[1]), get_n_eff(mt.pheno_data, x[0]))),
            hl.enumerate(mt.summary_stats),
        )
    )

    mt = mt.filter_entries(
        hl.is_defined(mt.summary_stats[P_FIELDS[0]]) &
        hl.is_defined(mt.summary_stats[P_FIELDS[1]]) &
        hl.is_defined(mt.summary_stats[P_FIELDS[2]])) # TODO: check if this is necessary?
    mt = mt.annotate_entries(**{p_field: hl.map(lambda x: _edit_pvalue(x), mt.summary_stats[p_field]) for p_field in P_FIELDS})

    for test in list(P_TESTS.keys()):
        print(f'Meta analyzing {test} results...')
        two_tail = test == 'Burden'
        beta_field = f'BETA_{test}'
        p_field = P_TESTS[test]

        if two_tail:
            mt = mt.annotate_entries(**{p_field: mt[p_field]/2},
                                     **{beta_field: mt.summary_stats[beta_field]})
        else:
            mt = mt.annotate_entries(**{beta_field: hl.map(lambda x: hl.int(hl.is_defined(x)), mt[p_field])})
        mt = mt.annotate_entries(**{f'weighted_Z_numerator_{test}' :
                                    hl.map(lambda x,y, z: hl.sqrt(x)*(-hl.qnorm(y))*  hl.sign(z),
                                           mt.summary_stats.N, mt[p_field], mt[beta_field])})
        mt = mt.annotate_entries(**{f'META_Stats_{test}': hl.sum(mt[f'weighted_Z_numerator_{test}']) / (hl.sqrt(hl.sum(mt.summary_stats.N)))}, )

        if two_tail:
            mt = mt.annotate_entries(**{f'META_Pvalue_{test}': 2*hl.pnorm(hl.abs(mt[f'META_Stats_{test}']), lower_tail=False)})
        else:
            mt = mt.annotate_entries(**{f'META_Pvalue_{test}': hl.pnorm(mt[f'META_Stats_{test}'], lower_tail=False)})

    return mt
