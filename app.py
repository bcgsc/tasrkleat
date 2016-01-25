import os
import re

import ruffus as R
import logging.config
logging.config.fileConfig('logging.config')
logger = logging.getLogger(__name__)

import utils as U
from config import CONFIG

import pprint
logger.info(pprint.pformat(CONFIG))

@R.mkdir(CONFIG['input_gs_bam'], R.formatter(),
         os.path.join(CONFIG['output_dir'], 'download_bam'))
@R.originate(
    os.path.join(CONFIG['output_dir'], 'download_bam', os.path.basename(CONFIG['input_gs_bam'])),
    # use extra param to store the flag filename
    [os.path.join(CONFIG['output_dir'], 'download_bam', 'download_bam.log'),
     os.path.join(CONFIG['output_dir'], 'download_bam', 'download_bam.SUCCESS')],
)
@U.timeit
def download_bam(output_bam, extras):
    log, flag = extras
    cmd = ('gsutil cp {bam} {outdir} 2>&1 '
           '| tee {log}').format(bam=CONFIG['input_gs_bam'],
                                 outdir=os.path.dirname(output_bam),
                                 log=log)
    U.execute(cmd, flag)


@R.mkdir(download_bam, R.formatter(), '{subpath[0][1]}/sort_bam_by_name')
@R.transform(download_bam, R.formatter(), [
    # name a is arbitarily selected, a short prefix for clarity and simplicity
    # since the sample can be disambiguated by the result path already. More
    # names will be selected in alphabetical order in the rest of the code.
    '{subpath[0][1]}/sort_bam_by_name/a.bam',
    '{subpath[0][1]}/sort_bam_by_name/sort_bam_by_name.log',
    '{subpath[0][1]}/sort_bam_by_name/sort_bam_by_name.SUCCESS'
])
@U.timeit
def sort_bam_by_name(input_bam, outputs):
    logger.info(outputs)
    bam, log, flag = outputs
    output_prefix = re.sub('\.bam$', '', bam)
    num_cpus = CONFIG['num_cpus']
    cmd = ('samtools sort -@ {num_cpus} -n {input_bam} {output_prefix} 2>&1 '
           '| tee {log}').format(**locals())
    U.execute(cmd, flag)


# @R.mkdir(sort_bam_by_name, R.formatter(), '{subpath[0][1]}/bam2fastq')
# @R.transform(sort_bam_by_name, R.formatter(), [
#     '{subpath[0][1]}/bam2fastq/a_1.fq',
#     '{subpath[0][1]}/bam2fastq/a_2.fq',
#     '{subpath[0][1]}/bam2fastq/bam2fastq.log',
#     '{subpath[0][1]}/bam2fastq/bam2fastq.SUCCESS'
# ])
# @U.timeit
# def bam2fastq(inputs, outputs):
#     input_bam, _, _ = inputs
#     fq1, fq2, log, flag = outputs
#     cmd = ('bedtools bamtofastq -i {input_bam} -fq {fq1} -fq2 {fq2} 2>&1 '
#            '| tee {log}').format(**locals())
#     U.execute(cmd, flag)

# instead of converting bam to fastq, which takes a lot of time, just remove
# the @CO lines from the bam, which causes problematic parsing in
# biobloomcategorizer, the same problem as described in
# https://groups.google.com/forum/#!msg/abyss-users/uDoJjgPWeu4/fJ-SYGN8XLsJ
@R.mkdir(sort_bam_by_name, R.formatter(), '{subpath[0][1]}/remove_CO_header')
@R.transform(sort_bam_by_name, R.formatter(), [
    '{subpath[0][1]}/remove_CO_header/a.bam',
    '{subpath[0][1]}/remove_CO_header/remove_CO_header.log',
    '{subpath[0][1]}/remove_CO_header/remove_CO_header.SUCCESS'
])
@U.timeit
def remove_CO_header(inputs, outputs):
    input_bam, _, _ = inputs
    output_bam, log, flag = outputs
    num_cpus = CONFIG['num_cpus']
    cmd = ('samtools view -h -@ {num_cpus} {input_bam} '
           '| grep -v "^@CO" '
           '| samtools view -Sb -@ {num_cpus} - > {output_bam} 2>&1'
           '| tee {log}').format(**locals())
    U.execute(cmd, flag)


@R.follows(remove_CO_header)
@R.mkdir(remove_CO_header, R.formatter(), '{subpath[0][1]}/download_bf')
@R.originate(
    os.path.join(CONFIG['output_dir'], 'download_bf', os.path.basename(CONFIG['input_gs_bf'])),
    [os.path.join(CONFIG['output_dir'], 'download_bf', 'download_bf.log'),
     os.path.join(CONFIG['output_dir'], 'download_bf', 'download_bf.SUCCESS')],
)
@U.timeit
def download_bf(output_bf, extras):
    log, flag = extras
    cmd = ('gsutil cp {bf} {outdir} 2>&1 '
           '| tee {log}').format(bf=CONFIG['input_gs_bf'],
                               outdir=os.path.dirname(output_bf), log=log)
    U.execute(cmd, flag)


@R.mkdir(download_bf, R.formatter(), '{subpath[0][1]}/extract_bf')
@R.transform(download_bf, R.formatter(), [
    # name p is arbitarily selected, a short prefix for clarity and simplicity,
    # the sample can be disambiguated by the result path already
    '{subpath[0][1]}/extract_bf/b.bf',
    '{subpath[0][1]}/extract_bf/b.txt',
    '{subpath[0][1]}/extract_bf/extract_bf.SUCCESS'
])
@U.timeit
def extract_bf(input_tar_gz, outputs):
    bf, txt, flag = outputs
    tar_gz_prefix = re.sub('\.tar\.gz$', '', os.path.basename(input_tar_gz))
    outdir = os.path.dirname(bf)
    cmd = ('tar zxf {input_tar_gz} -C {outdir} '
           '&& mv -v {outdir}/{tar_gz_prefix}/{tar_gz_prefix}.bf {outdir}/b.bf '
           '&& mv -v {outdir}/{tar_gz_prefix}/{tar_gz_prefix}.txt {outdir}/b.txt '
           '&& rmdir -v {outdir}/{tar_gz_prefix}').format(**locals())
    CONFIG['input_bf'] = bf
    U.execute(cmd, flag)


@R.follows(extract_bf)
@R.mkdir(remove_CO_header, R.formatter(), '{subpath[0][1]}/biobloomcategorizer')
@R.transform(remove_CO_header, R.formatter(), [
    '{{subpath[0][1]}}/biobloomcategorizer/a_{0}.fq.gz'.format(CONFIG['steps']['biobloomcategorizer']['bf_name']),
    '{{subpath[0][1]}}/biobloomcategorizer/p_{0}.fq.gz'.format(CONFIG['steps']['biobloomcategorizer']['bf_name']),
    '{subpath[0][1]}/biobloomcategorizer/biobloomcategorizer.log',
    '{subpath[0][1]}/biobloomcategorizer/biobloomcategorizer.SUCCESS'
])
@U.timeit
def biobloomcategorizer(inputs, outputs):
    input_bam, _, _ = inputs
    output_fq1, output_fq2, log, flag = outputs
    output_prefix = os.path.join(os.path.dirname(flag), 'a')
    bf = CONFIG['input_bf']
    num_cpus = CONFIG['num_cpus']
    # considered using -d, but then the paired-end reads get interlaced into a
    # single file, which would become problematic when the paired-end read
    # names aren't distinguishable
    cmd = ("biobloomcategorizer -p {output_prefix} -e -i -f '{bf}' -t {num_cpus} "
           "--gz_output --fq {input_bam} 2>&1 | tee {log}".format(**locals()))
    U.execute(cmd, flag)




                                            
if __name__ == "__main__":
    R.pipeline_run()
