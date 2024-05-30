# https://github.com/Nealelab/ukb_common/blob/9ec920ce9e0c39c4b9baf896ae820c9e0d7f167a/utils/annotations.py
import hail as hl


PLOF_CSQS = ["transcript_ablation", "splice_acceptor_variant",
             "splice_donor_variant", "stop_gained", "frameshift_variant"]

MISSENSE_CSQS = ["stop_lost", "start_lost", "transcript_amplification",
                 "inframe_insertion", "inframe_deletion", "missense_variant"]

SYNONYMOUS_CSQS = ["stop_retained_variant", "synonymous_variant"]

OTHER_CSQS = ["mature_miRNA_variant", "5_prime_UTR_variant",
              "3_prime_UTR_variant", "non_coding_transcript_exon_variant", "intron_variant",
              "NMD_transcript_variant", "non_coding_transcript_variant", "upstream_gene_variant",
              "downstream_gene_variant", "TFBS_ablation", "TFBS_amplification", "TF_binding_site_variant",
              "regulatory_region_ablation", "regulatory_region_amplification", "feature_elongation",
              "regulatory_region_variant", "feature_truncation", "intergenic_variant"]

# Unused: protein_altering_variant, incomplete_terminal_codon_variant, coding_sequence_variant
# TODO: question, what to do with: "splice_region_variant"
# TODO: question, "missense-damaging" vs "damaging_missense"


def annotation_case_builder(worst_csq_by_gene_canonical_expr, use_loftee: bool = True, use_polyphen_and_sift: bool = False,
                            strict_definitions: bool = False):
    case = hl.case(missing_false=True)
    if use_loftee:
        case = (case
                .when(worst_csq_by_gene_canonical_expr.lof == 'HC', 'pLoF')
                .when(worst_csq_by_gene_canonical_expr.lof == 'LC', 'LC'))
    else:
        case = case.when(hl.set(PLOF_CSQS).contains(worst_csq_by_gene_canonical_expr.most_severe_consequence), 'pLoF')
    if use_polyphen_and_sift:
        case = (case
                .when(missense.contains(mt.vep.worst_csq_for_variant_canonical.most_severe_consequence) &
                      (mt.vep.worst_csq_for_variant_canonical.polyphen_prediction == "probably_damaging") &
                      (mt.vep.worst_csq_for_variant_canonical.sift_prediction == "deleterious"), "damaging_missense")
                .when(missense.contains(mt.vep.worst_csq_for_variant_canonical.most_severe_consequence), "other_missense"))
    else:
        if strict_definitions:
            case = case.when(worst_csq_by_gene_canonical_expr.most_severe_consequence == 'missense_variant', 'missense')
        else:
            case = case.when(hl.set(MISSENSE_CSQS).contains(worst_csq_by_gene_canonical_expr.most_severe_consequence), 'missense')
    if strict_definitions:
        case = case.when(worst_csq_by_gene_canonical_expr.most_severe_consequence == 'synonymous_variant', 'synonymous')
    else:
        case = case.when(hl.set(SYNONYMOUS_CSQS).contains(worst_csq_by_gene_canonical_expr.most_severe_consequence), 'synonymous')
    case = case.when(hl.set(OTHER_CSQS).contains(worst_csq_by_gene_canonical_expr.most_severe_consequence), 'non-coding')
    return case.or_missing()


def annotation_case_builder_ukb_legacy(worst_csq_by_gene_canonical_expr):
    return (hl.case(missing_false=True)
            .when(worst_csq_by_gene_canonical_expr.lof == 'HC', 'pLoF')
            .when(worst_csq_by_gene_canonical_expr.lof == 'LC', 'LC')
            .when((worst_csq_by_gene_canonical_expr.most_severe_consequence == 'missense_variant') |
                  (worst_csq_by_gene_canonical_expr.most_severe_consequence == 'inframe_insertion') |
                  (worst_csq_by_gene_canonical_expr.most_severe_consequence == 'inframe_deletion'), 'missense')
            # TODO: add stop/start lost
            .when(worst_csq_by_gene_canonical_expr.most_severe_consequence == 'synonymous_variant', 'synonymous')
            .or_missing())


def create_gene_map_ht(context_ht, indel_ht, check_gene_contigs=False, freq_field=None):
    from gnomad.utils.vep import process_consequences

    def format_vep_ht(ht, freq_field, check_gene_contigs):
        fields = ['variant_id', 'gene_id', 'gene_symbol', 'annotation']
        if freq_field is not None:
            ht = ht.annotate(_af=ht[freq_field])
            fields.append('_af')
        ht = process_consequences(ht)
        ht = ht.explode(ht.vep.worst_csq_by_gene_canonical)
        ht = ht.filter(ht.vep.worst_csq_by_gene_canonical.gene_id.startswith('ENSG'))
        ht = ht.annotate(
            variant_id=ht.locus.contig + ':' + hl.str(ht.locus.position) + ':' + ht.alleles[0] + ':' + ht.alleles[1],
            annotation=annotation_case_builder(ht.vep.worst_csq_by_gene_canonical))
        if check_gene_contigs:
            gene_contigs = ht.group_by(
                gene_id=ht.vep.worst_csq_by_gene_canonical.gene_id,
                gene_symbol=ht.vep.worst_csq_by_gene_canonical.gene_symbol,
            ).aggregate(
                contigs=hl.agg.collect_as_set(ht.locus.contig)
            )
            assert gene_contigs.all(hl.len(gene_contigs.contigs) == 1)
        ht = ht.annotate(gene_id=ht.vep.worst_csq_by_gene_canonical.gene_id,
                         gene_symbol=ht.vep.worst_csq_by_gene_canonical.gene_symbol,)

        ht = ht.select(**{field: ht[field] for field in fields})
        return ht

    context_ht = format_vep_ht(context_ht, freq_field, check_gene_contigs)
    indel_ht = format_vep_ht(indel_ht, freq_field, check_gene_contigs)
    ht = context_ht.union(indel_ht)
    collect_field = (ht.variant_id, ht._af) if freq_field is not None else ht.variant_id
    gene_map_ht = ht.group_by(
        gene_id=ht.gene_id,
        gene_symbol=ht.gene_symbol,
    ).partition_hint(100).aggregate(
        interval=hl.interval(
            start=hl.locus(hl.agg.take(ht.locus.contig, 1)[0], hl.agg.min(ht.locus.position), reference_genome='GRCh38'),
            end=hl.locus(hl.agg.take(ht.locus.contig, 1)[0], hl.agg.max(ht.locus.position), reference_genome='GRCh38')
        ),
        variants=hl.agg.group_by(ht.annotation, hl.agg.collect(collect_field)),
    )
    return gene_map_ht


def tmp_annotation_case_builder(worst_csq_by_gene_canonical_expr, use_loftee: bool = True, use_polyphen_and_sift: bool = False,
                            strict_definitions: bool = False):
    case = hl.case(missing_false=True)
    if use_loftee:
        case = (case
                .when(worst_csq_by_gene_canonical_expr.lof == 'HC', 'pLoF')
                .when(worst_csq_by_gene_canonical_expr.lof == 'LC', 'LC'))
    else:
        case = case.when(hl.set(PLOF_CSQS).contains(worst_csq_by_gene_canonical_expr), 'pLoF')
    if use_polyphen_and_sift:
        case = (case
                .when(missense.contains(mt.vep.worst_csq_for_variant_canonical.most_severe_consequence) &
                      (mt.vep.worst_csq_for_variant_canonical.polyphen_prediction == "probably_damaging") &
                      (mt.vep.worst_csq_for_variant_canonical.sift_prediction == "deleterious"), "damaging_missense")
                .when(missense.contains(mt.vep.worst_csq_for_variant_canonical.most_severe_consequence), "other_missense"))
    else:
        if strict_definitions:
            case = case.when(worst_csq_by_gene_canonical_expr.most_severe_consequence == 'missense_variant', 'missense')
        else:
            case = case.when(hl.set(MISSENSE_CSQS).contains(worst_csq_by_gene_canonical_expr), 'missense')
    if strict_definitions:
        case = case.when(worst_csq_by_gene_canonical_expr.most_severe_consequence == 'synonymous_variant', 'synonymous')
    else:
        case = case.when(hl.set(SYNONYMOUS_CSQS).contains(worst_csq_by_gene_canonical_expr), 'synonymous')
    case = case.when(hl.set(OTHER_CSQS).contains(worst_csq_by_gene_canonical_expr), 'non-coding')
    return case.or_missing()


def create_tmp_gene_map_ht(ht, check_gene_contigs=False, freq_field=None):
    if freq_field is not None:
        ht = ht.annotate(_af=freq_field)
    ht = ht.annotate(gene_id = ht['gene_id'],
                     gene_symbol = ht['gene_symbol'])
    ht = ht.filter(ht.gene_id.startswith('ENSG'))
    ht = ht.annotate(
        variant_id=ht.locus.contig + ':' + hl.str(ht.locus.position) + ':' + ht.alleles[0] + ':' + ht.alleles[1],
        annotation=tmp_annotation_case_builder(ht['consequence'][0], use_loftee=False)
    )
    if check_gene_contigs:
        gene_contigs = ht.group_by('gene_id', 'gene_symbol').aggregate(
            contigs=hl.agg.collect_as_set(ht.locus.contig)
        )
        assert gene_contigs.all(hl.len(gene_contigs.contigs) == 1)

    collect_field = (ht.variant_id, ht._af) if freq_field is not None else ht.variant_id
    gene_map_ht = ht.group_by('gene_id', 'gene_symbol').partition_hint(100).aggregate(
        interval=hl.interval(
            start=hl.locus(hl.agg.take(ht.locus.contig, 1)[0], hl.agg.min(ht.locus.position)),
            end=hl.locus(hl.agg.take(ht.locus.contig, 1)[0], hl.agg.max(ht.locus.position))
        ),
        variants=hl.agg.group_by(ht.annotation, hl.array(hl.agg.collect_as_set(collect_field))),
    )
    return gene_map_ht


def post_process_gene_map_ht(gene_ht, freq_cutoff):
    print(f'Frequency cutoff: {freq_cutoff}')
    groups = ['pLoF', 'missense|LC', 'synonymous']
    variant_groups = hl.map(lambda group: group.split('\\|').flatmap(lambda csq: gene_ht.variants.get(csq)), groups)
    gene_ht = gene_ht.transmute(
        variant_groups=hl.zip(groups, variant_groups)
    ).explode('variant_groups')
    gene_ht = gene_ht.transmute(annotation=gene_ht.variant_groups[0],
                                variants=hl.sorted(gene_ht.variant_groups[1]))

    common_variants: hl.expr.ArrayExpression = gene_ht.variants.filter(lambda x: x[1] >= freq_cutoff)
    rare_variants = gene_ht.variants.filter(lambda x: x[1] < freq_cutoff)
    variants = common_variants.map(lambda x: (gene_ht.annotation, True, [x])).append((gene_ht.annotation, False, rare_variants))
    gene_ht = gene_ht.select('interval', variants=variants).explode('variants')
    gene_ht = gene_ht.transmute(annotation=gene_ht.variants[0], common_variant=gene_ht.variants[1], variants=hl.sorted(gene_ht.variants[2].map(lambda x: x[0])))
    gene_ht = gene_ht.filter(~gene_ht.common_variant)
    gene_ht = gene_ht.annotate(annotation = hl.if_else(gene_ht.annotation.matches('\\|'), gene_ht.annotation.replace('\\|', ''), gene_ht.annotation))
    gene_ht = gene_ht.annotate(annotation = (gene_ht.annotation+' ')*hl.len(gene_ht.variants))
    gene_ht = gene_ht.key_by(start=gene_ht.interval.start)
    return gene_ht.filter(hl.len(gene_ht.variants) > 0)