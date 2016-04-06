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
    os.path.join(CONFIG['output_dir'], 'abyss', 'cba-6.fa')
])
@U.timeit
def abyss(inputs, outputs):
    # better use absolute path because of -C option in abyss, which changes the
    # relative path.
    input_fq1, input_fq2 = [os.path.abspath(_) for _ in inputs]
    contigs_fa = outputs[0]
    cfg = CONFIG['steps']['abyss']
    too_small, read_count = U.fastq_too_small(input_fq1)

    if too_small:
        cfg.update(locals())
        msg = ('Only {read_count} (expect > {num_reads_cutoff}) reads are found '
               'in\n\t{input_fq1}\n\t{input_fq2}\n'
               'too small for assembly').format(**cfg)
        logging.info(msg)
        return

    outdir = os.path.dirname(contigs_fa)
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

# will replace abyss with the following transabyss commands
# transabyss --pe experiment/tasrkleat-results/biobloomcategorizer/cba_targetUTRcell2009_1.fq experiment/tasrkleat-results/biobloomcategorizer/cba_targetUTRcell2009_2.fq  --outdir k25  --name aaa --island 0 -c 1 --kmer 25
# transabyss --pe experiment/tasrkleat-results/biobloomcategorizer/cba_targetUTRcell2009_1.fq experiment/tasrkleat-results/biobloomcategorizer/cba_targetUTRcell2009_2.fq  --outdir k35  --name aaa --island 0 -c 1 --kmer 35
# transabyss --pe experiment/tasrkleat-results/biobloomcategorizer/cba_targetUTRcell2009_1.fq experiment/tasrkleat-results/biobloomcategorizer/cba_targetUTRcell2009_2.fq  --outdir k45  --name aaa --island 0 -c 1 --kmer 45
# # output ./transabyss-merged.fa
# transabyss-merge k{25,35,45}/aaa-6.fa --mink 25 --maxk 45


@R.mkdir(abyss, R.formatter(), '{subpath[0][1]}/align_contigs_to_genome')
@R.transform(abyss, R.formatter(), [
    '{subpath[0][1]}/align_contigs_to_genome/cba.sort.bam',
])
def align_contigs_to_genome(inputs, outputs):
    contigs_fa = inputs[0]
    output_bam = outputs[0]
    cfg = CONFIG['steps']['align_contigs_to_genome']
    cfg.update(locals())
    cmd = ('bwa mem {reference_genome_bwa_index} {contigs_fa} '
           '| samtools view -h -F 2052 -S - '
           '| samtools sort -o {output_bam} - '.format(**cfg))
    U.execute(cmd)


# left as a reference of converting bam to fa
# @R.mkdir(align_contigs_to_genome, R.formatter(), '{subpath[0][1]}/contigs_bam2fa')
# @R.transform(align_contigs_to_genome, R.formatter(), [
#     '{subpath[0][1]}/contigs_bam2fa/cba.fa',
# ])
# def contigs_bam2fa(inputs, outputs):
#     """converts aligned contigs to fasta format"""
#     input_bam = inputs[0]
#     output_fa = outputs[0]
#     cmd = 'samtools bam2fq {input_bam} | seqtk seq -A > {output_fa}'.format(**locals())
#     U.execute(cmd)


@R.transform(abyss, R.formatter(), [
    '{subpath[0][1]}/contigs_bam2fa/cba-6.fa.pac',
    '{subpath[0][1]}/contigs_bam2fa/cba-6.fa.bwt',
    '{subpath[0][1]}/contigs_bam2fa/cba-6.fa.ann',
    '{subpath[0][1]}/contigs_bam2fa/cba-6.fa.amb',
    '{subpath[0][1]}/contigs_bam2fa/cba-6.fa.sa',
])
def index_contigs_fa(inputs, outputs):
    contigs_fa = inputs[0]
    cmd = 'bwa index {contigs_fa}'.format(**locals())
    U.execute(cmd)


@R.follows(index_contigs_fa)
@R.mkdir(biobloomcategorizer, R.formatter(), '{subpath[0][1]}/align_reads_to_contigs')
@R.transform(biobloomcategorizer, R.formatter(), [
    '{subpath[0][1]}/align_reads_to_contigs/cba.bam',
])
def align_reads_to_contigs(inputs, outputs):
    input_fq1, input_fq2 = inputs
    # the same as the output contigs.fa from abyss
    index = os.path.join(CONFIG['output_dir'], 'abyss', 'cba-6.fa')
    output_bam = outputs[0]
    cmd = ('bwa mem {index} {input_fq1} {input_fq2} '
           '| samtools view -h -F 2052 -S - '
           '| samtools sort -o {output_bam}'.format(**locals()))
    U.execute(cmd)


if __name__ == "__main__":
    R.pipeline_run()
