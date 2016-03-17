#!/bin/env python

import ruffus as R
import multiprocessing


def parse_args():
    parser = R.cmdline.get_argparse(
        description='variant-call',
        usage='require python-2.7.x')

    parser.add_argument(
        '-i', '--input-bam', required=True,
        help='input bam file')

    parser.add_argument(
        '--input-ref-fa', required=True,
        help='reference genome in fasta format (e.g. hg19.fa)')

    parser.add_argument(
        '--input-ref-gff', required=True,
        help='reference gff file')

    parser.add_argument(
        '-t', '--num-threads', default=multiprocessing.cpu_count(),
        help=('if not specified, will use the number of available cpus on '
              'the machine'))

    parser.add_argument(
        '--output-log',
        help='output log file, default to <current_dir>/bamqc.log')
    
    args = parser.parse_args()
    return args
