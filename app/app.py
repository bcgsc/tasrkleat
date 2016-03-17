#!/usr/bin/env python

import os
import itertools

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


@R.transform(CONFIG['input_bam'], R.formatter(),
    '{path[0]}/{basename[0]}.parallel.sh'
)
def gen_parallel_cmd(input_bam, output_sh):
    output_dir = os.path.join(CONFIG['output_dir'], 'mpileup_vcf')
    cfg = CONFIG['steps']['gen_parallel_cmd']
    cfg.update(locals())
    cmd = ('vcfutils.pl '
           'splitchr '
           ' -l 50000000 '
           '{ref_fa}.fai '
           '| xargs '
           '-i echo "samtools mpileup -C50 -ABuf {ref_fa} -r {{}} {input_bam} | bcftools view -cvg - > {output_dir}/{{}}.vcf" '
           '> {output_sh}'.format(**cfg))
    U.execute(cmd)


@R.follows(index_bam)
@R.mkdir(CONFIG['input_bam'], R.formatter(),
         os.path.join(CONFIG['output_dir'], 'mpileup_vcf'))
@R.split(gen_parallel_cmd,
    os.path.join(CONFIG['output_dir'], 'mpileup_vcf', '*.vcf'),
)
@U.timeit
def mpileup_vcf(input_sh, output_vcfs):
    num_cpus = CONFIG['num_cpus']
    cmd = 'parallel -j {num_cpus} < {input_sh}'.format(**locals())
    U.execute(cmd)


@R.mkdir(mpileup_vcf, R.formatter(), '{subpath[0][1]}/find_hexamer_snvs')
@R.transform(mpileup_vcf, R.formatter(), [
    '{subpath[0][1]}/find_hexamer_snvs/{basename[0]}.tsv',
    # '{subpath[0][1]}/find_hexamer_snvs/{basename[0]}.sorted_filtered.vcf'
])
@U.timeit
def find_hexamer_snvs(input_vcf, outputs):
    output_tsv, _ = outputs
    cfg = CONFIG['steps']['find_hexamer_snvs']
    cfg.update(locals())
    cmd = ('find_hexamer_snvs.py '
           '{input_vcf} '
           '{ref_fa} '
           '{ref_gff} '
           '{output_tsv}'.format(**cfg))
    U.execute(cmd)


def sort_key(val):
    # e.g. val /experiment/mpileup-hexmer-snv-results/find_hexamer_snvs/9:100000001-141213431.tsv
    base = os.path.basename(val).split(':')[0]
    try:
        res = int(base)
    except ValueError:
        res = base
    return res


@R.mkdir(find_hexamer_snvs, R.formatter(),
    os.path.join(CONFIG['output_dir'], 'merge_tsv'))
@R.merge(find_hexamer_snvs,
    os.path.join(CONFIG['output_dir'], 'merge_tsv', 'merged.tsv'),
)
@U.timeit
def collate(input_tsvs, output_tsv):
    with open(output_tsv, 'wt') as opf:
        sorted_tsvs = sorted(itertools.chain(*input_tsvs), key=sort_key)
        for k, tsv in enumerate(sorted_tsvs):
            logger.info('working on {0}'.format(tsv))
            with open(tsv) as inf:
                if k != 0:
                    # skip header
                    inf.readline()
                for line in inf:
                    opf.write(line)


if __name__ == "__main__":
    R.pipeline_printout_graph(os.path.join(CONFIG['output_dir'], 'flowchart.svg'))
    R.pipeline_run(logger=logger, multiprocess=CONFIG['num_cpus'])
