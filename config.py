import os
import re
import multiprocessing

from argsparser import parse_args


def gen_config():
    args = parse_args()
    config = {
        'input_gs_bam': args.input_gs_bam,
        'input_gs_bf': args.input_gs_bf,
        'num_cpus': multiprocessing.cpu_count(),

        'steps': {
            'biobloomcategorizer': {
                'bf_name': re.sub('\.tar\.gz$', '', os.path.basename(args.input_gs_bf))
            },

            'abyss': {
                'num_reads_cutoff': args.abyss_num_reads_cutoff,
                'kmer_size': args.abyss_kmer_size,
            },

            'upload': {
                'output_gs_bucket': args.upload_gs_bucket,
            },
        },

    }


    if 'prefix' not in config or not config['prefix']:
        config['prefix'] = os.path.basename(config['input_gs_bam']).rstrip('.bam')

    if 'output_dir' not in config or not config['output_dir']:
        config['output_dir'] = os.path.join(os.getcwd(), 'tasrcloud_results', config['prefix'])

    try:
        os.makedirs(config['output_dir'])
    except OSError:
        pass

    config['logging'] = {
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
                 'filename': os.path.join(config['output_dir'], 'tasrcloud.log'),
                 'formatter': 'standard'
             },

         },

     }
    return config

CONFIG = gen_config()
