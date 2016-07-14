#!/usr/bin/env python

import os
import sys
import re
import itertools
import subprocess

import ruffus as R

from config import CONFIG
import logging.config
logging.config.dictConfig(CONFIG['logging'])
logger = logging.getLogger('tasrkleat')

import utils as U

import pprint
logger.info('\n{0}'.format(pprint.pformat(CONFIG)))


@R.mkdir(CONFIG['input_tar'], R.formatter(),
         os.path.join(CONFIG['output_dir'], 'extract_tarball'))
@R.split(CONFIG['input_tar'],
         os.path.join(CONFIG['output_dir'], 'extract_tarball', 'cba_[12].fastq'))
def extract_tarball(input_tarball, outputs):
    output_dir = os.path.join(CONFIG['output_dir'], 'extract_tarball')
    if input_tarball.endswith('.tar.gz'):
        cmd = 'tar zxfv {input_tarball} -C {output_dir}'.format(**locals())
    elif input_tarball.endswith('.tar'):
        cmd = 'tar  xfv {input_tarball} -C {output_dir}'.format(**locals())
    else:
        raise ValueError(
            'unrecognized tarball format "{0}"'.format(input_tarball))
    U.execute(cmd)

    # list of extracted files
    exfiles = [os.path.join(output_dir, __) for __ in os.listdir(output_dir)]
    in_fq1s = sorted([__ for __ in exfiles
                      if re.search('.*_1\.fastq(?:\.gz)?$', __)])
    in_fq2s = sorted([__ for __ in exfiles
                      if re.search('.*_2\.fastq(?:\.gz)?$', __)])

    out_fq1 = os.path.join(output_dir, 'cba_1.fastq')
    out_fq2 = os.path.join(output_dir, 'cba_2.fastq')
    if len(in_fq1s) == len(in_fq1s):
        if in_fq1s[0].endswith('.gz'):
            U.execute('gunzip -c {0} > {1}'.format(' '.join(in_fq1s), out_fq1))
            U.execute('gunzip -c {0} > {1}'.format(' '.join(in_fq2s), out_fq2))
        else:
            if len(in_fq1s) == 1:
                os.rename(in_fq1s[0], out_fq1)
                os.rename(in_fq2s[0], out_fq2)
            elif len(in_fq1s) > 1:
                U.execute('cat {0} > {1}'.format(' '.join(in_fq1s), out_fq1))
                U.execute('cat {0} > {1}'.format(' '.join(in_fq2s), out_fq2))
    else:
        raise ValueError("the numbers of _1.fastq and _2.fastq don't match: "
                         "{0} vs {1}".format(len(in_fq1s), len(in_fq2s)))


@R.mkdir(extract_tarball, R.formatter(), '{subpath[0][1]}/biobloomcategorizer')
@R.collate(extract_tarball, R.formatter(), [
    os.path.join(CONFIG['output_dir'],
                 'biobloomcategorizer',
                 'cba_{bf_name}_{i}.{ext}'.format(
                     bf_name=CONFIG['steps']['biobloomcategorizer']['bf_name'],
                     i=i,
                     ext=ext))
    for ext in ['fq', 'fq.gz']
    for i in [1, 2]
])
@U.timeit
def biobloomcategorizer(inputs, outputs):
    """categorize input reads with biobloomcategorizer"""
    input_fq1, input_fq2 = sorted(inputs)
    output_fq1, output_fq2, output_fq1_gz, output_fq2_gz = outputs
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

    # if put gzip part in a new task, it's too cumbersome to compose the output
    # file names as seen above
    U.execute('gzip -c {output_fq1} > {output_fq1_gz}'.format(**locals()))
    U.execute('gzip -c {output_fq2} > {output_fq2_gz}'.format(**locals()))


@R.mkdir(biobloomcategorizer, R.formatter(), '{subpath[0][1]}/transabyss')
@R.transform(biobloomcategorizer, R.formatter(), [
    os.path.join(CONFIG['output_dir'], 'transabyss', 'merged.fa')
])
@U.timeit
def transabyss(inputs, outputs):
    """Run transabyss for three kmer sizes and transabyss-merge"""
    # better use absolute path because of -C option in abyss, which changes the
    # relative path.
    input_fq1, input_fq2, _, _ = [os.path.abspath(_) for _ in inputs]
    contigs_fa = outputs[0]
    cfg = CONFIG['steps']['transabyss']
    kmer_sizes = cfg['kmer_sizes']
    too_small, read_count = U.fastq_too_small(input_fq1)

    cfg.update(locals())
    if too_small:
        msg = ('Only {read_count} (expect > {num_reads_cutoff}) reads are found '
               'in\n\t{input_fq1}\n\t{input_fq2}\n'
               'too small for assembly').format(**cfg)
        logging.warning(msg)

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

    kmer_sizes_to_merge = []
    fas_to_merge = []
    for ksize in kmer_sizes:
        fa = os.path.join(outdir, 'k{0}/aaa-6.fa'.format(ksize))
        if os.path.exists(fa):
            fas_to_merge.append(fa)
            kmer_sizes_to_merge.append(ksize)
    fas_to_merge = ' '.join(fas_to_merge)

    if len(kmer_sizes_to_merge) == 0:
        logger.error('no contig files (i.e. aaa-6.fa) are assembled for '
                         'any of the kmer sizes: {0} from {1} & {2}'.format(
                             ', '.join(map(str, kmer_sizes)),
                             input_fq1,
                             input_fq2))
        sys.exit(1)

    merge_cmd = ('transabyss-merge {fas_to_merge} '
                 '--mink {min_k} '
                 '--maxk {max_k} '
                 '--out {contigs_fa}'.format(
                     min_k=min(kmer_sizes_to_merge),
                     max_k=max(kmer_sizes_to_merge),
                     **locals()))
    U.execute(merge_cmd)


@R.mkdir(transabyss, R.formatter(), '{subpath[0][1]}/align_contigs2genome')
@R.transform(transabyss, R.formatter(), [
    '{subpath[0][1]}/align_contigs2genome/cba.sorted.bam',
])
def align_contigs2genome(inputs, outputs):
    contigs_fa = inputs[0]
    output_bam = outputs[0]
    output_bam_prefix = re.sub('\.bam$', '', output_bam)
    num_cpus = CONFIG['num_cpus']
    cfg = CONFIG['steps']['align_contigs2genome']
    cfg.update(locals())
    cmd = ('gmap '
           '--db hg19 '
           '--dir {reference_genome_gmap_index} '
           '--nthreads {num_cpus} '
           '--format samse '
           '--npaths 0 '
           '--chimera-margin 10 '
           '{contigs_fa}'
           '| samtools view -bhS -F 2052 - '
           # the api for samtools-0.1 and samtools-1.x are different!
           '| samtools sort - {output_bam_prefix}'.format(**cfg))
    U.execute(cmd)


@R.transform(align_contigs2genome, R.suffix('.bam'), output='.bam.bai')
def index_contigs2genome(inputs, output):
    input_bam = inputs[0]
    cmd = 'samtools index {input_bam}'.format(**locals())
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


@R.transform(transabyss, R.formatter(), [
    '{subpath[0][1]}/transabyss/merged.fa.pac',
    '{subpath[0][1]}/transabyss/merged.fa.bwt',
    '{subpath[0][1]}/transabyss/merged.fa.ann',
    '{subpath[0][1]}/transabyss/merged.fa.amb',
    '{subpath[0][1]}/transabyss/merged.fa.sa',
])
def index_contigs_fa(inputs, outputs):
    contigs_fa = inputs[0]
    cmd = 'bwa index {contigs_fa}'.format(**locals())
    U.execute(cmd)


@R.follows(index_contigs_fa)
@R.mkdir(biobloomcategorizer, R.formatter(), '{subpath[0][1]}/align_reads2contigs')
@R.transform(biobloomcategorizer, R.formatter(), [
    '{subpath[0][1]}/align_reads2contigs/cba.sorted.bam',
])
def align_reads2contigs(inputs, outputs):
    input_fq1, input_fq2, _, _ = inputs
    # the same as the output contigs.fa from abyss
    index = os.path.join(CONFIG['output_dir'], 'transabyss', 'merged.fa')
    output_bam = outputs[0]
    output_bam_prefix = re.sub('\.bam$', '', output_bam)
    cmd = ('bwa mem {index} {input_fq1} {input_fq2} '
           '| samtools view -bhS -F 2052 - '
           '| samtools sort - {output_bam_prefix}'.format(**locals()))
    U.execute(cmd)


@R.transform(align_reads2contigs, R.suffix('.bam'), output='.bam.bai')
def index_reads2contigs(inputs, output):
    input_bam = inputs[0]
    cmd = 'samtools index {input_bam}'.format(**locals())
    U.execute(cmd)


@R.mkdir(align_reads2contigs, R.formatter(), '{subpath[0][1]}/kleat')
@R.follows(index_reads2contigs)
@R.follows(index_contigs2genome)
@R.merge([transabyss,
          align_contigs2genome,
          align_reads2contigs],
         os.path.join(CONFIG['output_dir'], 'kleat/cba.KLEAT'))
def kleat(inputs, output):
    print(inputs)
    contigs_fa, c2g_bam, r2c_bam = list(itertools.chain(*inputs))
    output_prefix = re.sub('\.KLEAT$', '', output)
    cfg = CONFIG['steps']['kleat']
    cmd = ' '.join(['python',
                    # use the KLEAT.py that comes with package
                    os.path.join(os.path.dirname(__file__), 'KLEAT.py'),
                    c2g_bam,
                    contigs_fa,
                    cfg['reference_genome'],
                    cfg['annotations'],
                    r2c_bam,
                    output_prefix])
    U.execute(cmd)


def cleanup(outdir):
    """
    remove unwanted files before the results would be transferred to google
    cloud storage

    :outdir: the dir where initial tasrkleat results are located
    """
    # remove the original sequence files
    U.execute('rm -rfv {0}'.format(
        os.path.join(outdir, 'extract_tarball')))
    # remove unecessary fq files
    U.execute('rm -rfv {0}/*.fq'.format(
        os.path.join(outdir, 'biobloomcategorizer')))
    # remove huge hidden files in kleat. Use .[a-z][A-Z] instead of .* because
    # otherwise it will return non-zero for being unable to delete . and ..
    U.execute('rm -rfv {0}/.[a-zA-Z]*'.format(
        os.path.join(outdir, 'kleat')))


@R.follows(kleat)
def transfer():
    cleanup(CONFIG['output_dir'])

    cfg = CONFIG['steps']['transfer']
    cfg['output_dir'] = CONFIG['output_dir']
    cmd = 'gsutil cp -r {output_dir} {output_gsc_path}'.format(**cfg)
    # use subprocess.call instead of U.execute because the race between
    # tasrkleat writing gsutil's output to tasrkleat.log and gsutil uploading
    # tasrkleat.log to GCS can cause upload failure
    subprocess.call(cmd.split())


if __name__ == "__main__":
    # R.pipeline_printout_graph(
    #     os.path.join(CONFIG['output_dir'], 'flowchart.svg'), 'svg')
    R.cmdline.run(CONFIG['args'])
