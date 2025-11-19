#!/usr/bin/env python3

__author__ = "Wenhan Lu"

import hail as hl
import hailtop.batch as hb
import argparse



def main(args):
    hl.init(default_reference="GRCh38")

    if args.docker_image is None:
        image = hb.docker.build_python_image('us-central1-docker.pkg.dev/aou-neale-gwas/hail-batch/hail_python_3.9',
                                                        requirements=['hail'], python_version='3.9')
    else:
        image = args.docker_image

    backend = hb.ServiceBackend(
        billing_project=args.billing_project,
        remote_tmpdir=args.tmp_bucket
    )
    b = hb.Batch(
        name=args.batch_name,
        requester_pays_project=args.request_pays_project,
        default_python_image=image,
        backend=backend,
    )

    j = b.new_python_job(name=f"")
    j.call()
    # if args.test:
    #     break


    b.run()






if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--billing-project",
        help="Name of the billing project",
        nargs="?",
        default="all-by-aou"
    )
    parser.add_argument(
        "--tmp-bucket",
        help="Path to the temporary bucket",
        nargs="?",
        default='gs://aou_wlu/tmp'
    )
    parser.add_argument(
        "--batch-name",
        help="Name of the batch",
        nargs="?",
        default=''
    )
    parser.add_argument(
        "--request-pays-project",
        help="Name of the batch",
        nargs="?",
        default='aou-gwas-neale'
    )
    parser.add_argument(
        "--docker-image",
        nargs="?",
        help="Docker image to use",
        default='us-central1-docker.pkg.dev/aou-neale-gwas/hail-batch/wenhan-hail:0.2.117',
    )
    parser.add_argument(
        "--input-file-path",
        help="Path to the input file",
        nargs="?",
        default=""
    )
    parser.add_argument(
        "--output-file-path",
        help="Path to the output file",
        nargs="?",
        default=""
    )
    parser.add_argument(
        "--test",
        help="Whether to run test on a subset of the data",
        action="store_true"
    )
    parser.add_argument(
        "--update-table",
        help="Update some table",
        action="store_true"
    )
    parser.add_argument(
        "--data-type",
        nargs="?",
        help="Type of data being used",
        default="gene",
        choices=[""],
    )
    parser.add_argument('--apply_filter', help='Impose some filter', default=None, type=float)
    parser.add_argument('--overwrite', help='Whether to overwrite the existing file', action='store_true')
    args = parser.parse_args()

    main(args)


# hailctl auth login
# hailctl config set batch/tmp_dir gs://aou_wlu/tmp
# hailctl config set batch/remote_tmpdir gs://aou_wlu/tmp
# hailctl config set batch/billing_project all-by-aou
# hailctl config set query/backend batch (or spark)