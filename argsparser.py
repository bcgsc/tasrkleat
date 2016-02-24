import ruffus as R


def parse_args():
    parser = R.cmdline.get_argparse(
        description='tasrcloud',
        usage='require python-2.7.x')

    parser.add_argument(
        '-i', '--input-gs-bam', required=True,
        help=('input bam file from Google Cloud Storage (GCS), '
              'e.g. gs://<bucket-name>/<some.bam>'))

    parser.add_argument(
        '-o', '--upload-gs-bucket', required=True,
        help=('the GCS bucket to write back results'))

    parser.add_argument(
        '-r', '--refresh-token', required=True,
        help=('the path to a file containing the refresh token, which authorizes '
              'the code to read and write to relevant GCS bucket typically '
              'you\'re supposed to mount it somewhere in the file system when '
              'launch a job/jobs with Kubernetes'))

    parser.add_argument(
        '--output-log',
        help='output log file, default to <current_dir>/bamqc.log')
    
    args = parser.parse_args()
    return args
