import ruffus as R


def parse_args():
    parser = R.cmdline.get_argparse(
        description="tasrcloud",
        usage='require python-2.7.x')

    parser.add_argument(
        '-i', '--input-gs-bam', required=True,
        help=("input bam file from Google Cloud Storage (GCS), "
              "e.g. gs://<bucket-name>/<some.bam>"))
    
    parser.add_argument(
        '-f', '--input-gs-bf', required=True,
        help=("input bloomfilter from GCS, the bloomfilter (.bf) and its "
              "corresponding txt file (.txt) should be tar and gziped into "
              "a single .tar.gz file, e.g. gs://<bucket-name>/<some-name.tar.gz>, "
              "the name of the bloomfilter is assumed to be 'some-name'"))

    parser.add_argument(
        '-n', '--abyss-num-reads-cutoff', type=int, default=50,
        help=("Under which tasrcloud consider there are too few reads and "
              "run abyss at all. The default 50 is arbitrarily chosen"))

    parser.add_argument(
        '-k', '--abyss-kmer-size', type=int, required=True
        help=('the kmer size for running abyss, typically go for roughly '
              '~half the read length as a first trial. '
              'e.g. for 76bp, k=36 could be adequate'))

    parser.add_argument(
        '-o', '--upload-gs-bucket', required=True,
        help=("the GCS bucket to write back results"))

    parser.add_argument(
        '-r', '--refresh-token', required=True,
        help=("the path to a file containing the refresh token, which authorizes "
              "the code to read and write to relevant GCS bucket typically "
              "you're supposed to mount it somewhere in the file system when "
              "launch a job/jobs with Kubernetes"))
    
    args = parser.parse_args()
    return args
