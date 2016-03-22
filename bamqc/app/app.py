#!/usr/bin/env python

import os

import ruffus as R

from config import CONFIG
import logging.config
logging.config.dictConfig(CONFIG['logging'])
logger = logging.getLogger(__name__)

import utils as U

import pprint
logger.info('\n{0}'.format(pprint.pformat(CONFIG)))


@R.mkdir(CONFIG['input_bam'], R.formatter(),
         os.path.join(CONFIG['output_dir'], 'fastqc'))
@R.transform(CONFIG['input_bam'], R.formatter(), [
    os.path.join(CONFIG['output_dir'], 'fastqc/cba.html'),
    os.path.join(CONFIG['output_dir'], 'fastqc/cba.zip'),
])
@U.timeit
def fastqc(input_bam, outputs):
    output_html, output_zip = outputs
    output_dir = os.path.dirname(output_html)
    num_cpus = CONFIG['num_cpus']
    cmd = ("fastqc -o {output_dir} -t {num_cpus} {input_bam} "
           "&& mv -v {output_dir}/*.html {output_html} "
           "&& mv -v {output_dir}/*.zip  {output_zip} ".format(**locals()))
    U.execute(cmd)


if __name__ == "__main__":
    R.pipeline_run()
