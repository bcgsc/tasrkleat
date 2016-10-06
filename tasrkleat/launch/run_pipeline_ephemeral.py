#!/usr/bin/python

import sys
import os
import re
import math
import pprint
import time

from oauth2client.client import GoogleCredentials
from apiclient.discovery import build
import googleapiclient


def gen_id_path(input_file):
    return re.sub(r'^gs\:\/\/(?P<bucket_name>[^/]+)\/', '', input_file)


def gen_input_parameters():
    res = []

    for (name, path) in [
        ('input.tar', 'input.tar'),
        ('targets.bf', 'targets.bf'),
        ('targets.txt', 'targets.txt'),
        ('reference.fa', 'reference.fa'),
        ('reference.fa.fai', 'reference.fa.fai'),
        ('reference.gtf.gz', 'reference.gtf.gz'),
        ('reference.gtf.gz.tbi', 'reference.gtf.gz.tbi'),

        ('hg19.salcpchilddc', 'gmapdb/hg19.salcpchilddc'),
        ('hg19.salcpguide1024', 'gmapdb/hg19.salcpguide1024'),
        ('hg19.saindex64meta', 'gmapdb/hg19.saindex64meta'),
        ('hg19.sachildexc', 'gmapdb/hg19.sachildexc'),
        ('hg19.ref153offsets64strm', 'gmapdb/hg19.ref153offsets64strm'),
        ('hg19.ref153offsets64meta', 'gmapdb/hg19.ref153offsets64meta'),
        ('hg19.salcpexc', 'gmapdb/hg19.salcpexc'),
        ('hg19.genomebits128', 'gmapdb/hg19.genomebits128'),
        ('hg19.chromosome.iit', 'gmapdb/hg19.chromosome.iit'),
        ('hg19.contig', 'gmapdb/hg19.contig'),
        ('hg19.sarray', 'gmapdb/hg19.sarray'),
        ('hg19.chrsubset', 'gmapdb/hg19.chrsubset'),
        ('hg19.saindex64strm', 'gmapdb/hg19.saindex64strm'),
        ('hg19.chromosome', 'gmapdb/hg19.chromosome'),
        ('hg19.version', 'gmapdb/hg19.version'),
        ('hg19.ref153positions', 'gmapdb/hg19.ref153positions'),
        # empty dir
        # ('hg19.maps', 'gmapdb/hg19.maps'),
        ('hg19.sachildguide1024', 'gmapdb/hg19.sachildguide1024'),
        ('hg19.genomecomp', 'gmapdb/hg19.genomecomp'),
        ('hg19.contig.iit', 'gmapdb/hg19.contig.iit'),
    ]:
        res.append({
            'name': name,
            'localCopy': {
                'path': path,
                'disk': 'datadisk'
            }
        })

    for name in [
        'TRANSABYSS_KMER_SIZES', 'OUTPUT_GSC_PATH'
    ]:
        res.append({
            # localCopy unset, would be passed as-is as environmental
            # variable, so the name is capitalized
            'name': name,
        })
    return res


def gen_pipeline_args_inputs(
        input_gs_tar, transabyss_kmer_sizes, output_gsc_path):
    return {
        'input.tar': input_gs_tar,
        'targets.bf': 'gs://tasrkleat-static/targets.bf',
        'targets.txt': 'gs://tasrkleat-static/targets.txt',
        'reference.fa': 'gs://tasrkleat-static/hg19.fa',
        'reference.fa.fai': 'gs://tasrkleat-static/hg19.fa.fai',
        'reference.gtf.gz': 'gs://tasrkleat-static/ensembl.fixed.sorted.gz',
        'reference.gtf.gz.tbi': 'gs://tasrkleat-static/ensembl.fixed.sorted.gz.tbi',

        'hg19.chromosome': 'gs://tasrkleat-static/gmapdb/hg19.chromosome',
        'hg19.chromosome.iit': 'gs://tasrkleat-static/gmapdb/hg19.chromosome.iit',
        'hg19.chrsubset': 'gs://tasrkleat-static/gmapdb/hg19.chrsubset',
        'hg19.contig': 'gs://tasrkleat-static/gmapdb/hg19.contig',
        'hg19.contig.iit': 'gs://tasrkleat-static/gmapdb/hg19.contig.iit',
        'hg19.genomebits128': 'gs://tasrkleat-static/gmapdb/hg19.genomebits128',
        'hg19.genomecomp': 'gs://tasrkleat-static/gmapdb/hg19.genomecomp',
        'hg19.ref153offsets64meta': 'gs://tasrkleat-static/gmapdb/hg19.ref153offsets64meta',
        'hg19.ref153offsets64strm': 'gs://tasrkleat-static/gmapdb/hg19.ref153offsets64strm',
        'hg19.ref153positions': 'gs://tasrkleat-static/gmapdb/hg19.ref153positions',
        'hg19.sachildexc': 'gs://tasrkleat-static/gmapdb/hg19.sachildexc',
        'hg19.sachildguide1024': 'gs://tasrkleat-static/gmapdb/hg19.sachildguide1024',
        'hg19.saindex64meta': 'gs://tasrkleat-static/gmapdb/hg19.saindex64meta',
        'hg19.saindex64strm': 'gs://tasrkleat-static/gmapdb/hg19.saindex64strm',
        'hg19.salcpchilddc': 'gs://tasrkleat-static/gmapdb/hg19.salcpchilddc',
        'hg19.salcpexc': 'gs://tasrkleat-static/gmapdb/hg19.salcpexc',
        'hg19.salcpguide1024': 'gs://tasrkleat-static/gmapdb/hg19.salcpguide1024',
        'hg19.sarray': 'gs://tasrkleat-static/gmapdb/hg19.sarray',
        'hg19.version': 'gs://tasrkleat-static/gmapdb/hg19.version',

        'TRANSABYSS_KMER_SIZES': ' '.join(map(str, transabyss_kmer_sizes)),
        'OUTPUT_GSC_PATH': output_gsc_path,
    }


def create_pipeline_body(
        project_id,
        pipeline_name,

        input_gs_tar,
        output_bucket,
        transabyss_kmer_sizes,
        disk_size):
    """disk_size in GB"""
    id_path = gen_id_path(input_gs_tar)
    output_gsc_path = 'gs://{bucket}/{id_path}'.format(
        bucket=output_bucket, id_path=id_path)
    return {
        'ephemeralPipeline': {
            'projectId': project_id,
            'name': pipeline_name,
            'description': '',
            'docker': {
                'cmd': (
                    'app.py '
                    '--input-tar /mnt/data/input.tar '
                    '--input-bf /mnt/data/targets.bf '
                    '--transabyss-kmer-sizes ${TRANSABYSS_KMER_SIZES} '
                    '--reference-genome /mnt/data/reference.fa '
                    '--gtf /mnt/data/reference.gtf.gz '
                    '--reference-genome-gmap-index /mnt/data/gmapdb '
                    '--output-gsc-path ${OUTPUT_GSC_PATH}'
                ),

                # 'imageName': 'us.gcr.io/{0}/bamqc'.format(PROJECT_ID),
                'imageName': 'zyxue/test',
            },

            'inputParameters': gen_input_parameters(),

            'resources': {
                'disks': [
                    {
                        'name': 'datadisk',
                        'autoDelete': True,
                        'mountPoint': '/mnt/data',
                        'sizeGb': disk_size,
                        'type': 'PERSISTENT_HDD',
                    }
                ],

                # 'preemptible': True,
                'minimumCpuCores': 1,
                #  Default: 3.75 (GB)
                'minimumRamGb': 20,

                # 'zones': ["us-central1-a", "us-central1-b",
                #           "us-central1-c", "us-central1-f"]
            }
        },

        'pipelineArgs': {
            'inputs': gen_pipeline_args_inputs(
                input_gs_tar, transabyss_kmer_sizes, output_gsc_path),

            'logging': {
                'gcsPath': '{0}/logs'.format(output_gsc_path)
            },

            'projectId': PROJECT_ID,

            'serviceAccount': {
                'email': SERVICE_ACCOUNT,
                'scopes': [
                    'https://www.googleapis.com/auth/compute',
                    'https://www.googleapis.com/auth/devstorage.full_control',
                    'https://www.googleapis.com/auth/genomics',
                    # 'https://www.googleapis.com/auth/cloud-platform'
                ]
            }
        }
    }


if __name__ == '__main__':

    DEBUG = False

    # global variables below usually don't need change
    PIPELINE_NAME = 'tasrkleat'

    # OUTPUT_BUCKET = 'zx-googl-genomics-test'
    if DEBUG:
        OUTPUT_BUCKET = 'zx-trial'
    else:
        OUTPUT_BUCKET = 'tasrkleat'

    PROJECT_ID = 'isb-cgc-03-0006'
    SERVICE_ACCOUNT = '904618123662-compute@developer.gserviceaccount.com'

    credentials = GoogleCredentials.get_application_default()
    service = build('genomics', 'v1alpha2', credentials=credentials)

    # with open('bi_gcs_objects.csv') as inf:
        # for k, line in enumerate(inf):
        #     if k >= 50 and k < 200:
    # with open('bi_gcs_objects.fix.csv') as inf:
    #     for k, line in enumerate(inf):
    #         if True:
    # with open('gsc_gcs_objects.csv') as inf:
    #     for k, line in enumerate(inf):
    #         if k >= 700 and k < 1400:
    # with open('unc_gcs_objects.csv') as inf:
    #     for k, line in enumerate(inf):
    #         if k >= 7000 and k < 10000:

    with open('missing_runs_gcs_objects.csv') as inf:
        for k, line in enumerate(inf):
            if True:
                tar = line.split('\t')[0]
                size = math.ceil(float(line.split()[1]))
                length = int(line.split('\t')[2])

                # factor of 30 is experiential, additional 50GB for reference
                # data, etc.
                disk_size = size * 30 + 50

                if 45 <= length <= 50:
                    transabyss_kmer_sizes = [22, 32, 42]
                elif length in [75, 76]:
                    transabyss_kmer_sizes = [32, 52, 72]
                elif length == 100:
                    transabyss_kmer_sizes = [32, 62, 92]
                else:
                    raise ValueError('unsure of transabyss kmer sizes for read length: {0}'.format(length))

                print(tar, transabyss_kmer_sizes)

                body = create_pipeline_body(
                    PROJECT_ID, PIPELINE_NAME,
                    tar, OUTPUT_BUCKET,
                    transabyss_kmer_sizes,
                    disk_size)
                if DEBUG:
                    pprint.pprint(body)

                while True:
                    count = 0
                    try:
                        resp = service.pipelines().run(body=body).execute()
                        pprint.pprint(resp)
                        operation_id = resp['name']
                        with open('operations.txt', 'at') as opf:
                            rec = '{0}\t{1}\n'.format(tar, operation_id)
                            opf.write(rec)
                        break
                    except googleapiclient.errors.HttpError as err:
                        count += 1
                        print('{0}th retry due to {0}'.format(count, err))

                time.sleep(0.2)
