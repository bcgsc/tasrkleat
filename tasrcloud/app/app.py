#!/usr/bin/env python

import os
import re

import ruffus as R

from config import CONFIG
import logging.config
logging.config.dictConfig(CONFIG['logging'])
logger = logging.getLogger(__name__)

import utils as U

import pprint
logger.info('\n{0}'.format(pprint.pformat(CONFIG)))


@R.mkdir(CONFIG['input_fq'], R.formatter(),
         os.path.join(CONFIG['output_dir'], 'biobloomcategorizer'))
@R.collate(
    [CONFIG['input_fq'], CONFIG['input_fq2']],
    R.formatter(),
    [
        os.path.join(CONFIG['output_dir'], 'biobloomcategorizer',
                     'cba_{0}_1.fq'.format(CONFIG['steps']['biobloomcategorizer']['bf_name'])),
        os.path.join(CONFIG['output_dir'], 'biobloomcategorizer',
                     'cba_{0}_2.fq'.format(CONFIG['steps']['biobloomcategorizer']['bf_name']))
    ]
)
@U.timeit
def biobloomcategorizer(inputs, outputs):
    # input_fq1, input_fq2 = inputs
    # instead of the above line, assign input_fq1 & input_fq2, respectively
    input_fq1 = CONFIG['input_fq']
    input_fq2 = CONFIG['input_fq2']

    output_fq1, output_fq2 = outputs
    output_dir = os.path.dirname(output_fq1)
    output_prefix = os.path.join(output_dir, 'cba')
    num_cpus = CONFIG['num_cpus']
    cfg = CONFIG['steps']['biobloomcategorizer']
    cfg.update(locals())
    # considered using -d, but then the paired-end reads get interlaced into a
    # single file, which would become problematic when the paired-end read
    # names aren't distinguishable
    cmd = ("biobloomcategorizer "
           "-p {output_prefix} "
           "-e "
           "-i "
           "-f '{input_bf}' "
           "-t {num_cpus} "
           "--fq {input_fq1} {input_fq2}".format(**cfg))
    U.execute(cmd)


@R.mkdir(biobloomcategorizer, R.formatter(), '{subpath[0][1]}/abyss')
@R.transform(biobloomcategorizer, R.formatter(), [
    # '{subpath[0][1]}/abyss/coverage.hist',
    # '{subpath[0][1]}/abyss/cba-1.dot',
    # '{subpath[0][1]}/abyss/cba-1.fa',
    # '{subpath[0][1]}/abyss/cba-1.path',
    # '{subpath[0][1]}/abyss/cba-2.dot',
    # '{subpath[0][1]}/abyss/cba-2.dot1',
    # '{subpath[0][1]}/abyss/cba-2.fa',
    # '{subpath[0][1]}/abyss/cba-2.path',
    # '{subpath[0][1]}/abyss/cba-3.dist',
    # '{subpath[0][1]}/abyss/cba-3.dot',
    '{subpath[0][1]}/abyss/cba-3.fa',
    # '{subpath[0][1]}/abyss/cba-3.fa.fai',
    # '{subpath[0][1]}/abyss/cba-3.hist',
    # '{subpath[0][1]}/abyss/cba-4.dot',
    # '{subpath[0][1]}/abyss/cba-4.fa',
    # '{subpath[0][1]}/abyss/cba-4.fa.fai',
    # '{subpath[0][1]}/abyss/cba-4.path1',
    # '{subpath[0][1]}/abyss/cba-4.path2',
    # '{subpath[0][1]}/abyss/cba-4.path3',
    # '{subpath[0][1]}/abyss/cba-5.dot',
    # '{subpath[0][1]}/abyss/cba-5.fa',
    # '{subpath[0][1]}/abyss/cba-5.path',
    # '{subpath[0][1]}/abyss/cba-6.dist.dot',
    # '{subpath[0][1]}/abyss/cba-6.dot',
    '{subpath[0][1]}/abyss/cba-6.fa',
    # '{subpath[0][1]}/abyss/cba-6.hist',
    # '{subpath[0][1]}/abyss/cba-6.path',
    # '{subpath[0][1]}/abyss/cba-6.path.dot',
    # '{subpath[0][1]}/abyss/cba-7.dot',
    # '{subpath[0][1]}/abyss/cba-7.fa',
    # '{subpath[0][1]}/abyss/cba-7.path',
    # '{subpath[0][1]}/abyss/cba-8.dot',
    '{subpath[0][1]}/abyss/cba-8.fa',
    # '{subpath[0][1]}/abyss/cba-bubbles.fa',
    # # this is symlinked to cba-6.dot
    # '{subpath[0][1]}/abyss/cba-contigs.dot',
    # # this is symlinked to cba-6.fa
    # '{subpath[0][1]}/abyss/cba-contigs.fa',
    # '{subpath[0][1]}/abyss/cba-data',
    # '{subpath[0][1]}/abyss/cba-data.tar.gz',
    # '{subpath[0][1]}/abyss/cba-indel.fa',
    # # this is symlinked to cba-8.dot
    # '{subpath[0][1]}/abyss/cba-scaffolds.dot',
    # # this is symlinked to cba-8.fa
    # '{subpath[0][1]}/abyss/cba-scaffolds.fa',
    # # this is symlinked to cba-stats.tab
    # '{subpath[0][1]}/abyss/cba-stats',
    # '{subpath[0][1]}/abyss/cba-stats.csv',
    # '{subpath[0][1]}/abyss/cba-stats.md',
    # '{subpath[0][1]}/abyss/cba-stats.tab',
    # # this is symlinked to cba-3.fa
    # '{subpath[0][1]}/abyss/cba-unitigs.fa',
])
@U.timeit
def abyss(inputs, outputs):
    input_fq1, input_fq2 = inputs
    unitigs_fa, contigs_fa, scafflods_fa = outputs
    cfg = CONFIG['steps']['abyss']
    too_small, read_count = U.fastq_too_small(input_fq1)

    if too_small:
        cfg.update(locals())
        msg = ('Only {read_count} (expect > {num_reads_cutoff}) reads are found '
               'in\n\t{input_fq1}\n\t{input_fq2}\n'
               'too small for assembly').format(**cfg)
        logging.info(msg)
        with open(log, 'wt') as opf:
            opf.write('{0}\n'.format(msg))
        U.touch(flag)
        return

    outdir = os.path.dirname(unitigs_fa)
    num_cpus = CONFIG['num_cpus']

    # Note: name=a won't work for abyss-pe because of the particular way
    # how abyss reads command line parameters
    cmd = ("abyss-pe "
           "name=cba "
           "k={kmer_size} "
           "in='{input_fq1} {input_fq2}' "
           "np={num_cpus} "
           "-C {outdir}").format(kmer_size=cfg['kmer_size'], **locals())
    U.execute(cmd)


if __name__ == "__main__":
    R.pipeline_run()
