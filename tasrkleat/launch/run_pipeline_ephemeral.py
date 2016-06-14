#!/usr/bin/python

import sys
import os
import re
import pprint

from oauth2client.client import GoogleCredentials
from apiclient.discovery import build


PIPELINE_ID = '11869245215301103929' # 2 vCPU, 2GB

PROJECT_ID = 'isb-cgc-03-0006'
SERVICE_ACCOUNT = '904618123662-compute@developer.gserviceaccount.com'

# OUTPUT_BUCKET = 'zx-googl-genomics-test'
OUTPUT_BUCKET = 'zx-trial'

credentials = GoogleCredentials.get_application_default()
service = build('genomics', 'v1alpha2', credentials=credentials)


def gen_id_path(input_file):
    return re.sub(r'^gs\:\/\/(?P<bucket_name>[^/]+)\/', '', input_file)


def create_pipeline_body(input_gs_tar):
    id_path = gen_id_path(input_gs_tar)
    output_gsc_path = 'gs://{bucket}/{id_path}'.format(bucket=OUTPUT_BUCKET,
                                                       id_path=id_path)
    return {
        'ephemeralPipeline': {
            'projectId': PROJECT_ID,
            'name': 'tasrkleat',
            'description': '',
            'docker': {
                'cmd': (
                    # 'app.py '
                    # '--input-tar /mnt/data/input.tar '
                    # '--input-bf /mnt/data/input.bf '
                    # '--transabyss-kmer-sizes 22 32 '
                    # '--reference-genome experiment/hg19/hg19.fa '
                    # '--reference-genome-bwa-index experiment/hg19/bwa-index/hg19.fa '
                    # '--gtf experiment/KLEAT-2.5.0/ensembl.fixed.sorted.gz '

                    'app.py '
                    '--input-tar /mnt/data/input.tar '
                    '--input-bf /mnt/data/combined.bf '
                    '--transabyss-kmer-sizes ${TRANSABYSS_KMER_SIZES} '
                    '--reference-genome mnt/data/reference.fa '
                    '--reference-genome-bwa-index /mnt/data/reference.fa '
                    '--gtf reference.gtf.gz '
                ),

                # 'imageName': 'us.gcr.io/{0}/bamqc'.format(PROJECT_ID),
                'imageName': 'zyxue/test',
            },

            'inputParameters': [
                {
                    'name': 'inputTar',
                    'localCopy': {
                        'path': 'input.tar',
                        'disk': 'data'
                    }
                },

                {
                    'name': 'inputBf',
                    'localCopy': {
                        'path': 'combined.bf',
                        'disk': 'data'
                    }
                },

                {
                    'name': 'inputBfTxt',
                    'description': 'the .txt file that comes along with .bf',
                    'localCopy': {
                        'path': 'combined.txt',
                        'disk': 'data'
                    }
                },


                {
                    # localCopy unset, would be passed as-is as environmental
                    # variable, so the name is capitalized
                    'name': 'TRANSABYSS_KMER_SIZES',
                },

                {
                    'name': 'referenceGenome',
                    'localCopy': {
                        'path': 'reference.fa',
                        'disk': 'data'
                    }
                },

                {
                    'name': 'referenceGenomeFai',
                    'localCopy': {
                        'path': 'reference.fa.fai',
                        'disk': 'data'
                    }
                },

                # BEGIN bwa index files
                {
                    'name': 'ReferenceGenomeBwaIndexAnn',
                    'localCopy': {
                        'path': 'reference.fa.ann',
                        'disk': 'data'
                    }
                },


                {
                    'name': 'ReferenceGenomeBwaIndexBwt',
                    'localCopy': {
                        'path': 'reference.fa.bwt',
                        'disk': 'data'
                    }
                },


                {
                    'name': 'ReferenceGenomeBwaIndexPac',
                    'localCopy': {
                        'path': 'reference.fa.pac',
                        'disk': 'data'
                    }
                },


                {
                    'name': 'ReferenceGenomeBwaIndexSa',
                    'localCopy': {
                        'path': 'reference.fa.sa',
                        'disk': 'data'
                    }
                },

                {
                    'name': 'ReferenceGenomeBwaIndexAmb',
                    'localCopy': {
                        'path': 'reference.fa.amb',
                        'disk': 'data'
                    }
                },
                # END bwa index files 

                {
                    'name': 'referenceGtfGz',
                    'localCopy': {
                        'path': 'reference.gtf.gz',
                        'disk': 'data'
                    }
                },

                {
                    'name': 'referenceGtfGzTbi',
                    'localCopy': {
                        'path': 'reference.gtf.gz.tbi',
                        'disk': 'data'
                    }
                },

            ],

            'outputParameters': [
                {
                    'name': 'transabyssOutput',
                    'localCopy': {
                        'path': 'tasrkleat-results/transabyss/*',
                        'disk': 'data'
                    },
                },

                {
                    'name': 'kleatOutput',
                    'localCopy': {
                        'path': 'tasrkleat-results/kleat/*',
                        'disk': 'data'
                    },
                },

                {
                    'name': 'logOutput',
                    'localCopy': {
                        'path': 'tasrkleat-results/tasrkleat.log',
                        'disk': 'data'
                    },
                },

            ],

            'resources': {
                'disks': [
                    {
                        'name': 'data',
                        'autoDelete': True,
                        'mountPoint': '/mnt/data',
                        'sizeGb': 100,
                        'type': 'PERSISTENT_HDD',
                    }
                ],

                # 'preemptible': True,
                'minimumCpuCores': 1,
                # 'minimumRamGb': 25,
            }
        },

        'pipelineArgs': {
            'inputs': {
                'inputTar': input_gs_tar,
                'inputBf': 'gs://tasrkleat/static/combined.bf',
                'inputBfTxt': 'gs://tasrkleat/static/combined.txt',
                'TRANSABYSS_KMER_SIZES': '32',
                'referenceGenome': 'gs://tasrkleat/static/hg19.fa',
                'referenceGenomeFai': 'gs://tasrkleat/static/hg19.fa.fai',
                'ReferenceGenomeBwaIndexAnn': 'gs://tasrkleat/static/bwa-index/hg19.fa.ann',
                'ReferenceGenomeBwaIndexBwt': 'gs://tasrkleat/static/bwa-index/hg19.fa.bwt',
                'ReferenceGenomeBwaIndexPac': 'gs://tasrkleat/static/bwa-index/hg19.fa.pac',
                'ReferenceGenomeBwaIndexSa': 'gs://tasrkleat/static/bwa-index/hg19.fa.sa',
                'ReferenceGenomeBwaIndexAmb': 'gs://tasrkleat/static/bwa-index/hg19.fa.amb',
                'referenceGtfGz': 'gs://tasrkleat/static/ensembl.fixed.sorted.gz',
                'referenceGtfGzTbi': 'gs://tasrkleat/static/ensembl.fixed.sorted.gz.tbi',
            },

            'outputs': {
                'transabyssOutput': '{output_gsc_path}/transabyss/'.format(**locals()),
                'kleatOutput': '{output_gsc_path}/kleat/'.format(**locals()),
                'logOutput': '{output_gsc_path}/tasrkleat.log'.format(**locals()),
            },

            'logging': {
                'gcsPath': '{output_gsc_path}/logs'.format(**locals())
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
    tar = sys.argv[1]
    body = create_pipeline_body(tar)
    # pprint.pprint(body)
    resp = service.pipelines().run(body=body).execute()
    pprint.pprint(resp)
    operation_id = resp['name']
    with open('operations.txt', 'at') as opf:
        rec = '{0}\t{1}\n'.format(tar, operation_id)
        opf.write(rec)
