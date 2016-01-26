import os
import multiprocessing


CONFIG = {
    # must be a google cloud storage path
    'input_gs_bam': 'gs://tasrcloud-test-data/test.bam',
    # 'input_gs_bam': 'gs://isb-cgc-open-zxue/CESC/DNA-Seq/C836.HEC-108.1.bam',
    'input_gs_bf': 'gs://tasrcloud-bfs/targetUTRcell2009.tar.gz',
    'num_cpus': multiprocessing.cpu_count(),

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

try:
    os.makedirs(CONFIG['output_dir'])
except OSError:
    pass


CONFIG['logging'] = {
    'version': 1,
    'disable_existing_loggers': True,

     'loggers': {
         '__main__': {
             'handlers': ['screen', 'file'],
             'level': 'DEBUG',
             'propagate': True,
         },
         'utils': {
             'handlers': ['screen', 'file'],
             'level': 'DEBUG',
             'propagate': True,
         },
     },

     'formatters': {
         'verbose': {
             'format': '%(levelname)s|%(asctime)s|%(name)s|%(module)s|%(process)d|%(processName)s|%(relativeCreated)d|%(thread)d|%(threadName)s|%(msecs)d ms|%(pathname)s+%(lineno)d|%(funcName)s:%(message)s'
         },
         'standard': {
             'format': '%(levelname)s|%(asctime)s|%(name)s:%(message)s'
         },
         'colorful': {
             # https://github.com/borntyping/python-colorlog#with-dictconfig
             '()': 'colorlog.ColoredFormatter',
             'format': '%(log_color)s%(levelname)s%(reset)s|%(log_color)s[%(asctime)s]%(reset)s|%(log_color)s%(name)s%(reset)s:%(message)s'
         }
     },

     'handlers': {
         'screen':{
             'level': 'DEBUG',
             'class': 'logging.StreamHandler',
             'formatter': 'colorful'
         },
         'file': {
             'level': 'DEBUG',
             'class': 'logging.FileHandler',
             'filename': os.path.join(CONFIG['output_dir'], 'tasrcloud.log'),
             'formatter': 'standard'
         },

     },

 }
