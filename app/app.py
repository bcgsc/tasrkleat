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

@R.mkdir(CONFIG['input_gs_bam'], R.formatter(),
         os.path.join(CONFIG['output_dir'], 'download_bam'))
@R.originate(
    os.path.join(CONFIG['output_dir'], 'download_bam', os.path.basename(CONFIG['input_gs_bam'])),
    # use extra param to store the flag filename
    [os.path.join(CONFIG['output_dir'], 'download_bam', 'download_bam.log'),
     os.path.join(CONFIG['output_dir'], 'download_bam', 'download_bam.COMPLETE')],
)
@U.timeit
def download_bam(output_bam, extras):
    log, flag = extras
    cmd = ('{auth_gsutil} -m cp {bam} {outdir} 2>&1 | tee {log}').format(
        auth_gsutil=CONFIG['auth_gsutil'],
        bam=CONFIG['input_gs_bam'],
        outdir=os.path.dirname(output_bam),
        log=log)
    U.execute(cmd, flag)


# This task cannot be merged with download_bam because otherwise following
# tasks (e.g. bam2fastq will take gtf as a bam file
@R.mkdir(CONFIG['input_gs_gtf'], R.formatter(),
         os.path.join(CONFIG['output_dir'], 'download_gtf'))
@R.originate(
    os.path.join(CONFIG['output_dir'], 'download_gtf', os.path.basename(CONFIG['input_gs_gtf'])),
    [os.path.join(CONFIG['output_dir'], 'download_gtf', 'download_gtf.log'),
     os.path.join(CONFIG['output_dir'], 'download_gtf', 'download_gtf.COMPLETE')],
)
@U.timeit
def download_gtf(output_gtf, extras):
    log, flag = extras
    cmd = ('{auth_gsutil} -m cp {gtf} {outdir} 2>&1 | tee {log}').format(
        auth_gsutil=CONFIG['auth_gsutil'],
        gtf=CONFIG['input_gs_gtf'],
        outdir=os.path.dirname(output_gtf),
        log=log)
    U.execute(cmd, flag)


@R.mkdir(CONFIG['input_gs_star_index'], R.formatter(),
         os.path.join(CONFIG['output_dir'], 'download_star_index'))
@R.originate(
    os.path.join(CONFIG['output_dir'], 'download_star_index', os.path.basename(CONFIG['input_gs_star_index'])),
    [os.path.join(CONFIG['output_dir'], 'download_star_index', 'download_star_index.log'),
     os.path.join(CONFIG['output_dir'], 'download_star_index', 'download_star_index.COMPLETE')],
)
@U.timeit
def download_star_index(output_star_index, extras):
    log, flag = extras
    cmd = ('{auth_gsutil} -m cp -r {star_index} {outdir} 2>&1 | tee {log}').format(
        auth_gsutil=CONFIG['auth_gsutil'],
        star_index=CONFIG['input_gs_star_index'],
        outdir=os.path.dirname(output_star_index),
        log=log)
    U.execute(cmd, flag)


@R.mkdir(download_bam, R.formatter(), '{subpath[0][1]}/bam2fastq')
@R.transform(download_bam, R.formatter(), [
    '{subpath[0][1]}/bam2fastq/cba_1.fastq',
    '{subpath[0][1]}/bam2fastq/cba_2.fastq',
    '{subpath[0][1]}/bam2fastq/bam2fastq.log',
    '{subpath[0][1]}/bam2fastq/bam2fastq.COMPLETE'
])
@U.timeit
def bam2fastq(input_bam, outputs):
    fq1, fq2, log, flag = outputs
    output_dir = os.path.dirname(fq1)
    num_cpus = CONFIG['num_cpus']
    cmd = ('picard-tools SamToFastq I={input_bam} '
           'FASTQ={fq1} '
           'SECOND_END_FASTQ={fq2} '
           'UNPAIRED_FASTQ=discarded.fastq '
           '2>&1 | tee {log} '.format(**locals()))
    U.execute(cmd, flag)


@R.follows(download_gtf)
@R.follows(download_star_index)
@R.mkdir(bam2fastq, R.formatter(), '{subpath[0][1]}/star_align')
@R.transform(bam2fastq, R.formatter(), [
    '{subpath[0][1]}/star_align/cba.bam',
    # star align output example, for now just the bam and Log.out are needed
    # cba_Aligned.out.bam
    # cba_Log.final.out
    # cba_Log.out
    # cba_Log.progress.out
    # cba_SJ.out.tab
    # cba__STARgenome
    '{subpath[0][1]}/star_align/star_align.log',
    '{subpath[0][1]}/star_align/star_align.COMPLETE'
])
@U.timeit
def star_align(inputs, outputs):
    fq1, fq2, _, _ = inputs
    gtf = os.path.join(CONFIG['output_dir'], 'download_gtf',
                       os.path.basename(CONFIG['input_gs_gtf']))
    star_index = os.path.join(CONFIG['output_dir'], 'download_star_index',
                              os.path.basename(CONFIG['input_gs_star_index']))
    out_bam, log, flag = outputs
    output_dir = os.path.dirname(outputs[0])
    num_cpus = CONFIG['num_cpus']
    cmd = (
        'STAR '
        '--runThreadN {num_cpus} '
        '--sjdbGTFfile {gtf} '
        '--sjdbOverhang 100 '
        '--genomeDir {star_index}  '
        '--readFilesIn {fq1} {fq2} '
        '--alignIntronMin 30 '
        '--alignIntronMax 500000 '
        '--outFilterIntronMotifs RemoveNoncanonicalUnannotated '
        '--outFilterMismatchNmax 10 '
        '--outSAMstrandField intronMotif '
        '--genomeLoad NoSharedMemory '
        '--readMatesLengthsIn NotEqual '
        '--outSAMtype BAM Unsorted '
        '--outFileNamePrefix {output_dir}/cba_'
        '&& mv -v {output_dir}/cba_Log.out {log} 2>&1 | tee -a {output_dir}/cba_Log.out '
        '&& mv -v {output_dir}/cba_Aligned.out.bam {out_bam} 2>&1 | tee -a {log} '
        .format(**locals()))
    U.execute(cmd, flag)



@R.mkdir(star_align, R.formatter(), '{subpath[0][1]}/sort_bam')
@R.transform(star_align, R.formatter(), [
    '{subpath[0][1]}/sort_bam/cba.bam',
    '{subpath[0][1]}/sort_bam/sort_bam.log',
    '{subpath[0][1]}/sort_bam/sort_bam.COMPLETE'
])
@U.timeit
def sort_bam(inputs, outputs):
    input_bam, _, _ = inputs
    output_bam, log, flag = outputs
    output_bam_prefix = re.sub('\.bam$', '', output_bam)
    num_cpus = CONFIG['num_cpus']
    cmd = ('samtools sort -@ {num_cpus} {input_bam} {output_bam_prefix} '
           '2>&1 | tee {log} '.format(**locals()))
    U.execute(cmd, flag)


@R.mkdir(sort_bam, R.formatter(), '{subpath[0][1]}/upload')
@R.transform(sort_bam, R.formatter(), [
    '{subpath[0][1]}/upload/upload.log',
    '{subpath[0][1]}/upload/upload.COMPLETE',
])
@U.timeit
def upload(inputs, outputs):
    # e.g. /top_dir/first_subdir/second_subdir/somefile.log => /top_dir
    top_dir = re.search(r'/[^/]+', os.path.abspath(inputs[0])).group(0)
    # e.g. /top_dir => top_dir
    top_dir_name = top_dir.lstrip('/')
    log, flag = outputs
    cfg = CONFIG['steps']['upload']
    auth_gsutil=CONFIG['auth_gsutil']
    bucket_dir = os.path.join(cfg['output_gs_bucket'], top_dir_name)

    # Currently, there is no include, so has to exclude as many big files as
    # possible,
    # http://stackoverflow.com/questions/34111297/how-to-include-file-in-gsutil-rsync
    re_files_to_exclude = '|'.join([r'.*\.gtf$',
                                    r'.*\.fastq$',
                                    r'.*\.fq$',
                                    r'.*\.fasta$',
                                    r'.*\.fa$',
                                    r'.*\.fai$',
                                    r'.*\.dict$',
                                    ] + map(os.path.basename,
                                            [CONFIG['input_gs_bam']]) + [
            # ignore the whole directory
            r'.*{0}.*'.format(os.path.basename(CONFIG['input_gs_star_index']))
            ]
                                   )
    cmd = ("{auth_gsutil} -m rsync -c -x '{re_files_to_exclude}' -r "
           "{top_dir} {bucket_dir} "
           "2>&1 | tee {log}").format(**locals())
    U.execute(cmd, flag)


@R.follows(upload)
@U.timeit
def cleanup():
    """remove the output folder to save disk space"""
    cmd = "rm -rfv {0}".format(CONFIG['output_dir'])
    U.execute(cmd)


if __name__ == "__main__":
    R.pipeline_printout_graph("flowchart.svg", "svg", ['upload'],
                              user_colour_scheme = {"colour_scheme_index" :6})
    R.pipeline_run()
