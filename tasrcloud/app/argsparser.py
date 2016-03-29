import ruffus as R
import multiprocessing


def parse_args():
    parser = R.cmdline.get_argparse(
        description="tasrcloud",
        usage='require python-2.7.x')

    # parser.add_argument(
    #     '-i', '--input-bam', required=True,
    #     help=("input bam file"))


    parser.add_argument(
        '--input-fq', required=True,
        help=("input fastq file"))

    parser.add_argument(
        '--input-fq2', required=True,
        help=("input fastq file of the second end"))

    parser.add_argument(
        '-f', '--input-bf', required=True,
        help=("input bloomfilter: the bloomfilter (.bf) and its "
              "corresponding txt file (.txt)"))

    parser.add_argument(
        '-u', '--abyss-num-reads-cutoff', type=int, default=50,
        help=("Under which tasrcloud consider there are too few reads and "
              "run abyss at all. The default 50 is arbitrarily chosen"))

    parser.add_argument(
        '-k', '--abyss-kmer-size', type=int, required=True,
        help=('the kmer size for running abyss, typically go for roughly '
              '~half the read length as a first trial. '
              'e.g. for 76bp, k=36 could be adequate'))

    parser.add_argument(
        '-t', '--num-threads', default=multiprocessing.cpu_count(),
        help=('if not specified, will use the number of available cpus on '
              'the machine'))

    parser.add_argument(
        '--output-log',
        help='output log file, default to <current_dir>/bamqc.log')

    args = parser.parse_args()
    return args
