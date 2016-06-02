import ruffus as R
import multiprocessing


def parse_args():
    parser = R.cmdline.get_argparse(
        description="tasrkleat",
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
        '-u', '--num-reads-too-small-cutoff', type=int, default=50,
        help=("Under which it's considered that there are too few reads and "
              "run assembly at all. The default 50 is arbitrarily chosen"))

    parser.add_argument(
        '-k', '--transabyss-kmer-sizes', type=int, nargs=3, required=True,
        help=('the 3 kmer sizes for running transabyss, typically go for roughly '
              'e.g. --transabyss-kmer-sizes 25 35 45'))

    parser.add_argument(
        '--reference-genome', required=True,
        help=('should have been indexed already in the form of '
              'e.g. /path/to/hg19.fa.fai'))

    parser.add_argument(
        '--reference-genome-bwa-index', required=True,
        help=('e.g. /path/to/hg19/bwa-index/hg19.fa'))

    parser.add_argument(
        '--gtf', required=True,
        help=('genome annotations. e.g. /path/to/Homo_sapiens.GRCh37.75.gtf, '
              'and it must have been indexed'))

    # parser.add_argument(
    #     '-t', '--num-threads', default=multiprocessing.cpu_count(),
    #     help=('if not specified, will use the number of available cpus on '
    #           'the machine'))

    parser.add_argument(
        '--output-log',
        help='output log file, default to <current_dir>/bamqc.log')

    args = parser.parse_args()
    return args
