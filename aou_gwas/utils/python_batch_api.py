def load_jobs_by_batch_ids(batch_ids, billing_project: str = 'gnomad-production'):
    from hailtop.batch_client.client import BatchClient
    bc = BatchClient(billing_project=billing_project)
    if isinstance(batch_ids, int):
        batch_ids = [batch_ids]
    all_jobs = []
    for batch_id in batch_ids:
        batch = bc.get_batch(batch_id)
        all_jobs.extend(list(batch.jobs()))
    return all_jobs

jobs = load_jobs_by_batch_ids(['4612056', '4626052'], billing_project='gnomad-production')
# [..., {'batch_id': 4626052, 'job_id': 544, 'name': 'Run_VerifyBamID_HGDP01181.alt_bwamem_GRCh38DH.20181023.Yi', 'user': 'wlu', 'billing_project': 'gnomad-production', 'state': 'Success', 'exit_code': 0, 'duration': 2386174, 'cost': 0.18503390871418057, 'msec_mcpu': 0}, ...]
sum([job['duration'] for job in jobs if job['duration'] is not None])
# 2419865786
sum([job['cost'] for job in jobs if job['cost'] is not None])
# 155.2454743630818
sum([job['msec_mcpu'] for job in jobs if job['msec_mcpu'] is not None])