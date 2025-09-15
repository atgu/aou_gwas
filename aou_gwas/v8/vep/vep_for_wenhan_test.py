import hailtop.batch as hb
backend = hb.ServiceBackend()
b = hb.Batch(backend=backend)
j = b.new_bash_job()
j.image('us-central1-docker.pkg.dev/aou-neale-gwas/vep/vep_105_wenhan:v2.4')
j.command('true')
b.run()