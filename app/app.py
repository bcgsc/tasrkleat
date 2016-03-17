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


@R.transform(CONFIG['input_bam'], R.formatter(),
    '{path[0]}/{basename[0]}.bam.bai'
)
@U.timeit
def index_bam(input_bam, output_bai):
    cmd = 'samtools index {input_bam}'.format(**locals())
    U.execute(cmd)


@R.follows(index_bam)
@R.mkdir(CONFIG['input_bam'], R.formatter(), '{path[0]}/mpileup_vcf')
@R.transform(CONFIG['input_bam'], R.formatter(),
    '{path[0]}/mpileup_vcf/{basename[0]}.vcf'
)
@U.timeit
def mpileup_vcf(input_bam, output_vcf):
    cfg = CONFIG['steps']['mpileup_vcf']
    cfg.update(locals())
    cmd = ('samtools mpileup '
           '-C50 '
           '-ABuf '
           '{ref_fa} '
           '{input_bam} '
           '| bcftools view '
           '-cvg '
           '- > {output_vcf}'.format(**cfg))
    U.execute(cmd)


@R.mkdir(mpileup_vcf, R.formatter(), '{subpath[0][1]}/find_hexamer_snvs')
@R.transform(mpileup_vcf, R.formatter(),
    '{subpath[0][1]}/find_hexamer_snvs/{basename[0]}.tsv',
)
@U.timeit
def find_hexamer_snvs(input_vcf, output_tsv):
    cfg = CONFIG['steps']['find_hexamer_snvs']
    cfg.update(locals())
    cmd = ('find_hexamer_snvs.py '
           '{input_vcf} '
           '{ref_fa} '
           '{ref_gff} '
           '{output_tsv}'.format(**cfg))
    U.execute(cmd)


if __name__ == "__main__":
    R.pipeline_printout_graph(os.path.join(CONFIG['output_dir'], 'flowchart.svg'))
    R.pipeline_run(logger=logger)
