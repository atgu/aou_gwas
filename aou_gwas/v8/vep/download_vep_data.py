import hailtop.batch as hb


with hb.ServiceBackend('all-by-aou') as backend:
    # FIXME: switch 95 to 105    
    # https://ftp.ensembl.org/pub/release-105/variation/indexed_vep_cache/homo_sapiens_vep_105_GRCh38.tar.gz
    root = 'homo_sapiens_vep_105_GRCh38'
    tar_path = f'{root}.tar.gz'
    dest_bucket = 'vep_105'
    remote_vep_data_dest = f'gs://{dest_bucket}/'

    b = hb.Batch(backend=backend, name='download_vep_data')

    j1 = b.new_job()
    j1.image('google/cloud-sdk:slim') # before was ubuntu:20.04
    j1.storage('100Gi')
    j1.cpu(8)
    j1.command('apt-get update && apt-get install -y wget')
    
    # FIXME: switch 95 to 105
    j1.command(
        f'wget -q -c ftp://ftp.ensembl.org/pub/release-105/variation/indexed_vep_cache/{tar_path} -O - > {j1.vep_data_tar}')

    j1.command(f'mkdir -p {j1.vep_data}')
    j1.command(f'tar -xz -f {j1.vep_data_tar} -C {j1.vep_data}')

    j1.command(
        'for f in $(find ' + j1.vep_data + ' -name "*.csi"); do '
        '  rel="${f#' + j1.vep_data + '/}"; '
        '  gsutil cp "$f" "' + remote_vep_data_dest + '${rel}"; '
        'done'
    )

    j1.command(
        f'wget -q -c ftp://ftp.ensembl.org/pub/release-105/fasta/homo_sapiens/dna_index/Homo_sapiens.GRCh38.dna.toplevel.fa.gz -O - > {j1.vep_data}/homo_sapiens/105_GRCh38/Homo_sapiens.GRCh38.dna.toplevel.fa.gz')
    j1.command(
        f'wget -q -c ftp://ftp.ensembl.org/pub/release-105/fasta/homo_sapiens/dna_index/Homo_sapiens.GRCh38.dna.toplevel.fa.gz.fai -O - > {j1.vep_data}/homo_sapiens/105_GRCh38/Homo_sapiens.GRCh38.dna.toplevel.fa.gz.fai')
    j1.command(
        f'wget -q -c ftp://ftp.ensembl.org/pub/release-105/fasta/homo_sapiens/dna_index/Homo_sapiens.GRCh38.dna.toplevel.fa.gz.gzi -O - > {j1.vep_data}/homo_sapiens/105_GRCh38/Homo_sapiens.GRCh38.dna.toplevel.fa.gz.gzi')

    b.write_output(j1.vep_data, remote_vep_data_dest)
    b.run()


# gsutil -u aou-neale-gwas ls gs://hail-qob-vep-grch38-us-central1/
# gs://hail-qob-vep-grch38-us-central1/95_GRCh38_indexed.tar
# gs://hail-qob-vep-grch38-us-central1/gerp_conservation_scores.homo_sapiens.GRCh38.bw
# gs://hail-qob-vep-grch38-us-central1/human_ancestor.fa.gz
# gs://hail-qob-vep-grch38-us-central1/human_ancestor.fa.gz.fai
# gs://hail-qob-vep-grch38-us-central1/human_ancestor.fa.gz.gzi
# gs://hail-qob-vep-grch38-us-central1/loftee.sql
# gs://hail-qob-vep-grch38-us-central1/homo_sapiens/
# gs://hail-qob-vep-grch38-us-central1/homo_sapiens_backup/
# gsutil -u aou-neale-gwas cp gs://hail-qob-vep-grch38-us-central1/loftee.sql gs://vep_105/
# gsutil -u aou-neale-gwas cp gs://hail-qob-vep-grch38-us-central1/human_ancestor* gs://vep_105/
# gsutil -u aou-neale-gwas cp gs://hail-qob-vep-grch38-us-central1/gerp_conservation_scores.homo_sapiens.GRCh38.bw gs://vep_105/
