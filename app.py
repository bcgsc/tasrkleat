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

@R.follows(download_bam)
@R.mkdir(download_bam, R.formatter(), '{subpath[0][1]}/download_bf')
@R.originate(
    os.path.join(CONFIG['output_dir'], 'download_bf', os.path.basename(CONFIG['input_gs_bf'])),
    [os.path.join(CONFIG['output_dir'], 'download_bf', 'download_bf.log'),
     os.path.join(CONFIG['output_dir'], 'download_bf', 'download_bf.COMPLETE')],
)
@U.timeit
def download_bf(output_bf, extras):
    log, flag = extras
    cmd = ('gsutil cp {bf} {outdir} 2>&1 | tee {log}').format(
        bf=CONFIG['input_gs_bf'], outdir=os.path.dirname(output_bf), log=log)
    U.execute(cmd, flag)


@R.mkdir(download_bf, R.formatter(), '{subpath[0][1]}/extract_bf')
@R.transform(download_bf, R.formatter(), [
    '{subpath[0][1]}/extract_bf/cba.bf',
    '{subpath[0][1]}/extract_bf/cba.txt',
    '{subpath[0][1]}/extract_bf/extract_bf.COMPLETE'
])
@U.timeit
def extract_bf(input_tar_gz, outputs):
    bf, txt, flag = outputs
    tar_gz_prefix = re.sub('\.tar\.gz$', '', os.path.basename(input_tar_gz))
    outdir = os.path.dirname(bf)
    cmd = ('tar zxf {input_tar_gz} -C {outdir} '
           '&& mv -v {outdir}/{tar_gz_prefix}/{tar_gz_prefix}.bf {outdir}/cba.bf '
           '&& mv -v {outdir}/{tar_gz_prefix}/{tar_gz_prefix}.txt {outdir}/cba.txt '
           '&& rmdir -v {outdir}/{tar_gz_prefix}').format(**locals())
    CONFIG['input_bf'] = bf
    U.execute(cmd, flag)


@R.follows(extract_bf)
@R.mkdir(download_bam, R.formatter(), '{subpath[0][1]}/biobloomcategorizer')
@R.transform(download_bam, R.formatter(), [
    # because of .format(), {{}} means it's literal
    '{{subpath[0][1]}}/biobloomcategorizer/cba_{0}_1.fq'.format(CONFIG['steps']['biobloomcategorizer']['bf_name']),
    '{{subpath[0][1]}}/biobloomcategorizer/cba_{0}_2.fq'.format(CONFIG['steps']['biobloomcategorizer']['bf_name']),
    '{subpath[0][1]}/biobloomcategorizer/biobloomcategorizer.log',
    '{subpath[0][1]}/biobloomcategorizer/biobloomcategorizer.COMPLETE'
])
@U.timeit
def biobloomcategorizer(input_bam, outputs):
    output_fq1, output_fq2, log, flag = outputs
    output_dir = os.path.dirname(output_fq1)
    output_prefix = os.path.join(output_dir, 'cba')
    bf = CONFIG['input_bf']
    num_cpus = CONFIG['num_cpus']
    # considered using -d, but then the paired-end reads get interlaced into a
    # single file, which would become problematic when the paired-end read
    # names aren't distinguishable
    cmd = ("biobloomcategorizer -p {output_prefix} -e -i -f '{bf}' -t {num_cpus} "
           "--fq {input_bam} 2>&1 | tee {log}".format(**locals()))
    U.execute(cmd, flag)
    # don't delete since filesystem in the container is ephemeral anyway, code
    # left for ref
    # for f in os.listdir(output_dir):
    #     path_f = os.path.join(output_dir, f)
    #     if path_f not in outputs:
    #         logger.debug('removing {0}'.format(path_f))
    #         os.remove(path_f)

@R.mkdir(biobloomcategorizer, R.formatter(), '{subpath[0][1]}/abyss')
@R.transform(biobloomcategorizer, R.formatter(), [
    # '{subpath[0][1]}/abyss/coverage.hist',
    # '{subpath[0][1]}/abyss/cba-1.fa',
    # '{subpath[0][1]}/abyss/cba-bubbles.fa',
    # '{subpath[0][1]}/abyss/cba-1.adj',
    # '{subpath[0][1]}/abyss/cba-1.path',
    # '{subpath[0][1]}/abyss/cba-2.path',
    # '{subpath[0][1]}/abyss/cba-2.adj',
    # '{subpath[0][1]}/abyss/cba-3.adj',
    # '{subpath[0][1]}/abyss/cba-3.fa',
    # '{subpath[0][1]}/abyss/cba-indel.fa',
    # # this is a symlink to cba-3.fa
    # '{subpath[0][1]}/abyss/cba-unitigs.fa',
    # '{subpath[0][1]}/abyss/cba-stats.tab',
    # # this is a symlink to cba-stats.tab
    # '{subpath[0][1]}/abyss/cba-stats',
    # '{subpath[0][1]}/abyss/cba-stats.csv',
    # '{subpath[0][1]}/abyss/cba-stats.md',
    '{subpath[0][1]}/abyss/abyss.log',
    '{subpath[0][1]}/abyss/abyss.COMPLETE'
])
@U.timeit
def abyss(inputs, outputs):
    input_fq1, input_fq2, _, _ = inputs
    log, flag = outputs
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

    outdir = os.path.dirname(log)
    num_cpus = CONFIG['num_cpus']
    # as a note: name=a won't work for abyss-pe because of the particular way
    # how abyss reads command line parameters
    cmd = ("abyss-pe name=cba k={kmer_size} in='{input_fq1} {input_fq2}' "
           "np={num_cpus} 2>&1 -C {outdir} 2>&1 | tee {log}").format(
               kmer_size=cfg['kmer_size'], **locals())
    U.execute(cmd, flag)


@R.follows(abyss)
@R.mkdir(abyss, R.formatter(), '{subpath[0][1]}/upload')
@R.transform(abyss, R.formatter(), [
    '{subpath[0][1]}/upload/upload.log',
    '{subpath[0][1]}/upload/upload.COMPLETE',
])
@U.timeit
def upload(inputs, outputs):
    # e.g. /tasrcloud_results/sample_name/abyss/cba.log => /tasrcloud_results/sample_name
    input_dir = os.path.dirname(os.path.dirname(inputs[0]))
    log, flag = outputs
    top_log = CONFIG['logging']['handlers']['file']['filename']
    cfg = CONFIG['steps']['upload']
    auth_gsutil=CONFIG['auth_gsutil']

    # e.g. /tasrcloud_results/sample_name => sample_name
    sample_name = os.path.basename(input_dir)
    bucket_dir = os.path.join(cfg['output_gs_bucket'], sample_name)

    re_files_to_exclude = '|'.join(['.*\.fq$', '.*\.bam$'])
    cmd = ("{auth_gsutil} -m rsync -x '{re_files_to_exclude}' -r -d "
           "{input_dir} {bucket_dir}").format(**locals())
    U.execute(cmd, flag)


if __name__ == "__main__":
    R.pipeline_run()
