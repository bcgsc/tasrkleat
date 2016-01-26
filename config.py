import os
import multiprocessing


CONFIG = {
    # must be a google cloud storage path
    'input_gs_bam': 'gs://tasrcloud-test-data/test.bam',
    # 'input_gs_bam': 'gs://isb-cgc-open-zxue/CESC/DNA-Seq/C836.HEC-108.1.bam',
    'input_gs_bf': 'gs://tasrcloud-bfs/targetUTRcell2009.tar.gz',
    'num_cpus': multiprocessing.cpu_count(),
    
    'output_dir': '/tasrcloud_results',

    'steps': {
        'biobloomcategorizer': {
            'bf_name': 'targetUTRcell2009',
        },

        'abyss': {
            'kmer_size': 15,
        },
    },

}


if 'prefix' not in CONFIG or not CONFIG['prefix']:
    CONFIG['prefix'] = os.path.basename(CONFIG['input_gs_bam']).rstrip('.bam')

if 'output_dir' not in CONFIG or not CONFIG['output_dir']:
    CONFIG['output_dir'] = os.path.join(os.getcwd(), 'tasrcloud_results', CONFIG['prefix'])

