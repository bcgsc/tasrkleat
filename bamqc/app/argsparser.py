import ruffus as R
import multiprocessing


def parse_args():
    parser = R.cmdline.get_argparse(
        description='bamqc',
        usage='require python-2.7.x')

    parser.add_argument(
        '-i', '--input-bam', required=True,
        help='input bam file')

    parser.add_argument(
        '-t', '--num-threads', default=multiprocessing.cpu_count(),
        help=('if not specified, will use the number of available cpus on '
              'the machine'))

    parser.add_argument(
        '--output-log',
        help='output log file, default to <current_dir>/bamqc.log')
    
    args = parser.parse_args()
    return args
