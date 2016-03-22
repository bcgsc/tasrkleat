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


@R.mkdir(CONFIG['input_gs_ref_fa'], R.formatter(),
         os.path.join(CONFIG['output_dir'], 'download_ref_fa'))
@R.originate(
    os.path.join(CONFIG['output_dir'], 'download_ref_fa', os.path.basename(CONFIG['input_gs_ref_fa'])),
    [os.path.join(CONFIG['output_dir'], 'download_ref_fa', 'download_ref_fa.log'),
     os.path.join(CONFIG['output_dir'], 'download_ref_fa', 'download_ref_fa.COMPLETE')],
)
@U.timeit
def download_ref_fa(output_ref_fa, extras):
    log, flag = extras
    cmd = ('{auth_gsutil} -m cp {ref_fa} {outdir} 2>&1 | tee {log}').format(
        auth_gsutil=CONFIG['auth_gsutil'],
        ref_fa=CONFIG['input_gs_ref_fa'],
        outdir=os.path.dirname(output_ref_fa),
        log=log)
    U.execute(cmd, flag)


@R.transform(download_ref_fa, R.formatter(), [
    # fai has to be saved to the same directory for use
    '{subpath[0][1]}/download_ref_fa/{basename[0]}.fa.fai',
    '{subpath[0][1]}/download_ref_fa/index_ref_fa.log',
    '{subpath[0][1]}/download_ref_fa/index_ref_fa.COMPLETE'
])
@U.timeit
def index_ref_fa(input_fa, outputs):
    _, log, flag = outputs      # the name of fai is not necessary to be stored
    cmd = 'samtools faidx {input_fa} 2>&1 | tee {log}'.format(**locals())
    U.execute(cmd, flag)


# dict has to be in the same directory of the fa file, so it has to run after
# download_ref_fa
@R.follows(download_ref_fa)
@R.originate(
    os.path.join(CONFIG['output_dir'], 'download_ref_fa', os.path.basename(CONFIG['input_gs_ref_dict'])),
    [os.path.join(CONFIG['output_dir'], 'download_ref_fa', 'download_ref_dict.log'),
     os.path.join(CONFIG['output_dir'], 'download_ref_fa', 'download_ref_dict.COMPLETE')],
)
@U.timeit
def download_ref_dict(output_ref_dict, extras):
    log, flag = extras
    cmd = ('{auth_gsutil} -m cp {ref_dict} {outdir} 2>&1 | tee {log}').format(
        auth_gsutil=CONFIG['auth_gsutil'],
        ref_dict=CONFIG['input_gs_ref_dict'],
        outdir=os.path.dirname(output_ref_dict),
        log=log)
    U.execute(cmd, flag)


@R.mkdir(download_bam, R.formatter(), '{subpath[0][1]}/add_rg')
@R.transform(download_bam, R.formatter(), [
    '{subpath[0][1]}/add_rg/cba.bam',
    '{subpath[0][1]}/add_rg/add_rg.log',
    '{subpath[0][1]}/add_rg/add_rg.COMPLETE'
])
@U.timeit
def add_rg(input_bam, outputs):
    output_bam, log, flag = outputs
    output_bam_prefix = re.sub('\.bam$', '', output_bam)
    num_cpus = CONFIG['num_cpus']
    cmd = (
        'picard-tools AddOrReplaceReadGroups '
        'I={input_bam} '
        'O={output_bam} '
        'SO=coordinate '
        'RGID=id '
        'RGLB=library '
        'RGPL=platform '
        'RGPU=machine '
        'RGSM=sample '
        'TMP_DIR=/tmp '
        '2>&1 | tee {log}'.format(**locals()))
    U.execute(cmd, flag)


@R.mkdir(add_rg, R.formatter(), '{subpath[0][1]}/mark_dup')
@R.transform(add_rg, R.formatter(), [
    '{subpath[0][1]}/mark_dup/cba.bam',
    '{subpath[0][1]}/mark_dup/cba.metrics',
    '{subpath[0][1]}/mark_dup/mark_dup.log',
    '{subpath[0][1]}/mark_dup/mark_dup.COMPLETE'
])
@U.timeit
def mark_dup(inputs, outputs):
    input_bam, _, _ = inputs
    output_bam, output_metrics, log, flag = outputs
    output_bam_prefix = re.sub('\.bam$', '', output_bam)
    num_cpus = CONFIG['num_cpus']
    cmd = (
        'picard-tools MarkDuplicates '
        'I={input_bam} '
        'O={output_bam} '
        'CREATE_INDEX=true '
        'VALIDATION_STRINGENCY=SILENT '
        'M={output_metrics} '
        'TMP_DIR=/tmp '
        '2>&1 | tee {log}'.format(**locals()))
    U.execute(cmd, flag)


@R.follows(index_ref_fa)          # this guarantees @R.follows(download_ref_fa)
@R.follows(download_ref_dict)
@R.mkdir(mark_dup, R.formatter(), '{subpath[0][1]}/split_n_cigar_reads')
@R.transform(mark_dup, R.formatter(), [
    '{subpath[0][1]}/split_n_cigar_reads/cba.bam',
    '{subpath[0][1]}/split_n_cigar_reads/split_n_cigar_reads.log',
    '{subpath[0][1]}/split_n_cigar_reads/split_n_cigar_reads.COMPLETE'
])
@U.timeit
def split_n_cigar_reads(inputs, outputs):
    input_bam, _, _, _ = inputs
    output_bam, log, flag = outputs
    output_bam_prefix = re.sub('\.bam$', '', output_bam)
    num_cpus = CONFIG['num_cpus']
    ref_fa = os.path.join(CONFIG['output_dir'], 'download_ref_fa',
                          os.path.basename(CONFIG['input_gs_ref_fa']))
    cmd = (
        'java '
        '-jar /GenomeAnalysisTK.jar '
        '-T SplitNCigarReads '
        '-R {ref_fa} '
        '-I {input_bam} '
        '-o {output_bam} '
        '-rf ReassignOneMappingQuality '
        '-RMQF 255 '
        '-RMQT 60 '
        '-U ALLOW_N_CIGAR_READS '
        '2>&1 | tee {log} '.format(**locals()))
    U.execute(cmd, flag)


@R.mkdir(split_n_cigar_reads, R.formatter(), '{subpath[0][1]}/call_haplotype')
@R.transform(split_n_cigar_reads, R.formatter(), [
    '{subpath[0][1]}/call_haplotype/cba.vcf',
    '{subpath[0][1]}/call_haplotype/call_haplotype.log',
    '{subpath[0][1]}/call_haplotype/call_haplotype.COMPLETE'
])
@U.timeit
def call_haplotype(inputs, outputs):
    input_bam,  _, _ = inputs
    output_vcf, log, flag = outputs
    num_cpus = CONFIG['num_cpus']
    ref_fa = os.path.join(CONFIG['output_dir'], 'download_ref_fa',
                          os.path.basename(CONFIG['input_gs_ref_fa']))
    cmd = (
        'java '
        '-jar /GenomeAnalysisTK.jar '
        '-T HaplotypeCaller '
        '-nct {num_cpus} '
        '-R {ref_fa} '
        '-I {input_bam} '
        '-dontUseSoftClippedBases '
        '-stand_call_conf 20.0 '
        '-stand_emit_conf 20.0 '
        '-o {output_vcf} '
        '2>&1 | tee {log} '.format(**locals()))
    U.execute(cmd, flag)


@R.mkdir(call_haplotype, R.formatter(), '{subpath[0][1]}/filter_vcf')
@R.transform(call_haplotype, R.formatter(), [
    '{subpath[0][1]}/filter_vcf/cba.vcf',
    '{subpath[0][1]}/filter_vcf/filter_vcf.log',
    '{subpath[0][1]}/filter_vcf/filter_vcf.COMPLETE'
])
@U.timeit
def filter_vcf(inputs, outputs):
    input_vcf,  _, _ = inputs
    output_vcf, log, flag = outputs
    num_cpus = CONFIG['num_cpus']
    ref_fa = os.path.join(CONFIG['output_dir'], 'download_ref_fa',
                          os.path.basename(CONFIG['input_gs_ref_fa']))
    cmd = (
        'java '
        '-jar /GenomeAnalysisTK.jar '
        '-T VariantFiltration '
        '-R {ref_fa} '
        '-V {input_vcf} '
        '-window 35 '
        '-cluster 3 '
        '-filterName FS '
        '-filter "FS > 30.0" '
        '-filterName QD '
        '-filter "QD < 2.0" '
        '-o {output_vcf} '
        '2>&1 | tee {log} '.format(**locals()))
    U.execute(cmd, flag)


@R.mkdir(filter_vcf, R.formatter(), '{subpath[0][1]}/upload')
@R.transform(filter_vcf, R.formatter(), [
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
    re_files_to_exclude = '|'.join([r'.*\.bam$',
                                    r'.*\.gtf$',
                                    r'.*\.fastq$',
                                    r'.*\.fq$',
                                    r'.*\.fasta$',
                                    r'.*\.fa$',
                                    r'.*\.fai$',
                                    r'.*\.dict$',
                                    ] + map(
            os.path.basename, [CONFIG[_] for _ in ['input_gs_bam',
                                                   'input_gs_ref_fa',
                                                   'input_gs_ref_dict']]))
    cmd = ("{auth_gsutil} -m rsync -x '{re_files_to_exclude}' -r "
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
    R.pipeline_printout_graph("flowchart.svg", "svg", ['cleanup'],
                              user_colour_scheme = {"colour_scheme_index" :6})
    R.pipeline_run()

