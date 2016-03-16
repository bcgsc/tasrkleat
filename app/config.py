import os
import re
import multiprocessing

from argsparser import parse_args


def gen_config():
    project_id = os.environ.get('PROJECT_ID')
    if project_id is None:
        raise ValueError(
            'No PROJECT_ID environmental variable found, analysis aborted.')
        
    args = parse_args()

    config = {
        'input_bam': args.input_bam,

        # there is no need to put the input file inside the steps configuration
        # since it will be passed along by ruffus
        'steps': {
            'mpileup_vcf': {
                'ref_fa': args.input_ref_fa
            },

            'find_hexamer_snvs': {
                'ref_fa': args.input_ref_fa,
                'ref_gff': args.input_ref_gff
            },
        }

    }
        # 'input_gff': args.ref_gff,

    if args.num_threads:
        config['num_cpus'] = args.num_threads
    else:
        config['num_cpus'] = multiprocessing.cpu_count()

    # trying to mimic the directory hierarchy as the input by remove gs://
    # at the beginning and .bam at the end of the url
    config['prefix'] = re.sub('\.bam$', '', config['input_bam'])

    config['output_dir'] = os.path.join(
        # use - instead of _ to avoid confusion because when downloading
        # the zip from GCS, / will replaced with _
        os.getcwd(), config['prefix'], '{0}-results'.format(project_id))

    output_log_file = args.output_log
    if not output_log_file:
        output_log_file = os.path.join(
            config['output_dir'], '{0}.log'.format(project_id))

    try:
        os.makedirs(config['output_dir'])
    except OSError:
        pass

    config['logging'] = configure_logging_dict(output_log_file)
    return config


def configure_logging_dict(output_log_file):
    return {
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
                'filename': output_log_file,
                'formatter': 'standard'
                },
            },
        }

CONFIG = gen_config()
