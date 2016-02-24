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


@R.mkdir(download_bam, R.formatter(), '{subpath[0][1]}/fastqc')
@R.transform(download_bam, R.formatter(), [
    '{subpath[0][1]}/fastqc/cba.html',
    '{subpath[0][1]}/fastqc/fastqc.log',
    '{subpath[0][1]}/fastqc/fastqc.COMPLETE'
])
@U.timeit
def fastqc(input_bam, outputs):
    output_html, log, flag = outputs
    output_dir = os.path.dirname(output_html)
    num_cpus = CONFIG['num_cpus']
    cmd = ("fastqc -o {output_dir} -t {num_cpus} {input_bam} 2>&1 "
           "| tee {log} "
           "&& mv -v {output_dir}/*.html {output_dir}/cba.html "
           "| tee -a {log}".format(**locals()))
    U.execute(cmd, flag)

@R.mkdir(fastqc, R.formatter(), '{subpath[0][1]}/upload')
@R.transform(fastqc, R.formatter(), [
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

    # no need to upload the zip file from fastqc, and input bam
    re_files_to_exclude = '|'.join(['.*\.zip$', '.*\.bam$'])
    cmd = ("{auth_gsutil} -m rsync -x '{re_files_to_exclude}' -r -d "
           "{top_dir} {bucket_dir}").format(**locals())

    # -x '' would be invalid and cause error
    # cmd = ("{auth_gsutil} -m rsync -r -d "
    #        "{top_dir} {bucket_dir}").format(**locals())

    U.execute(cmd, flag)


if __name__ == "__main__":
    R.pipeline_run()
