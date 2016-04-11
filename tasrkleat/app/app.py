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


@R.mkdir(biobloomcategorizer, R.formatter(), '{subpath[0][1]}/transabyss')
@R.transform(biobloomcategorizer, R.formatter(), [
    os.path.join(CONFIG['output_dir'], 'transabyss', 'merged.fa')
])
@U.timeit
def transabyss(inputs, outputs):
    """Run transabyss for three kmer sizes and transabyss-merge"""
    # better use absolute path because of -C option in abyss, which changes the
    # relative path.
    input_fq1, input_fq2 = [os.path.abspath(_) for _ in inputs]
    contigs_fa = outputs[0]
    cfg = CONFIG['steps']['transabyss']
    kmer_sizes = cfg['kmer_sizes']
    too_small, read_count = U.fastq_too_small(input_fq1)

    cfg.update(locals())
    if too_small:
        msg = ('Only {read_count} (expect > {num_reads_cutoff}) reads are found '
               'in\n\t{input_fq1}\n\t{input_fq2}\n'
               'too small for assembly').format(**cfg)
        logging.info(msg)
        return

    outdir = os.path.dirname(contigs_fa)
    num_cpus = CONFIG['num_cpus']

    # Note: name=a won't work for abyss-pe because of the particular way
    # how abyss reads command line parameters
    for k in kmer_sizes:
        kmer_outdir = os.path.join(outdir, 'k{0}'.format(k))
        if not os.path.exists(kmer_outdir):
            os.mkdir(kmer_outdir)
        cmd = ('transabyss '
               '--pe {input_fq1} {input_fq2} '
               '--outdir {kmer_outdir} '
               '--name aaa '
               '--island 0 '
               '-c 1 --kmer {k}'.format(kmer_outdir=kmer_outdir, k=k, **cfg))
        U.execute(cmd)

    fas_to_merge = ' '.join([
        os.path.join(outdir, 'k{0}/aaa-6.fa'.format(__))
        for __ in kmer_sizes
    ])
    merge_cmd = ('transabyss-merge {fas_to_merge} '
                 '--mink {min_k} '
                 '--maxk {max_k} '
                 '--out {contigs_fa}'.format(
                     min_k=min(kmer_sizes),
                     max_k=max(kmer_sizes),
                     **locals()))
    U.execute(merge_cmd)


@R.mkdir(transabyss, R.formatter(), '{subpath[0][1]}/align_contigs_to_genome')
@R.transform(transabyss, R.formatter(), [
    '{subpath[0][1]}/align_contigs_to_genome/cba.sort.bam',
    '{subpath[0][1]}/align_contigs_to_genome/cba.sort.bam.bai',
])
def align_contigs_to_genome_and_index(inputs, outputs):
    contigs_fa = inputs[0]
    output_bam = outputs[0]
    cfg = CONFIG['steps']['align_contigs_to_genome']
    cfg.update(locals())
    cmd = ('bwa mem {reference_genome_bwa_index} {contigs_fa} '
           '| samtools view -h -F 2052 -S - '
           '| samtools sort -o {output_bam} - '.format(**cfg))
    U.execute(cmd)
    index_cmd = 'samtools index {output_bam}'.format(**cfg)
    U.execute(index_cmd)

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


@R.transform(transabyss, R.formatter(), [
    '{subpath[0][1]}/contigs_bam2fa/merged.fa.pac',
    '{subpath[0][1]}/contigs_bam2fa/merged.fa.bwt',
    '{subpath[0][1]}/contigs_bam2fa/merged.fa.ann',
    '{subpath[0][1]}/contigs_bam2fa/merged.fa.amb',
    '{subpath[0][1]}/contigs_bam2fa/merged.fa.sa',
])
def index_contigs_fa(inputs, outputs):
    contigs_fa = inputs[0]
    cmd = 'bwa index {contigs_fa}'.format(**locals())
    U.execute(cmd)


@R.follows(index_contigs_fa)
@R.mkdir(biobloomcategorizer, R.formatter(), '{subpath[0][1]}/align_reads_to_contigs')
@R.transform(biobloomcategorizer, R.formatter(), [
    '{subpath[0][1]}/align_reads_to_contigs/cba.bam',
    '{subpath[0][1]}/align_reads_to_contigs/cba.bam.bai',
])
def align_reads_to_contigs_and_index(inputs, outputs):
    input_fq1, input_fq2 = inputs
    # the same as the output contigs.fa from abyss
    index = os.path.join(CONFIG['output_dir'], 'abyss', 'cba-6.fa')
    output_bam = outputs[0]
    cmd = ('bwa mem {index} {input_fq1} {input_fq2} '
           '| samtools view -h -F 2052 -S - '
           '| samtools sort -o {output_bam}'.format(**locals()))
    U.execute(cmd)
    index_cmd = 'samtools index {output_bam}'.format(**locals())
    U.execute(index_cmd)


if __name__ == "__main__":
    R.pipeline_run()
