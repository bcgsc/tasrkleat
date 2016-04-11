import os
import re
import multiprocessing

from argsparser import parse_args


def gen_config():
    project_id = os.environ['PROJECT_ID']
    args = parse_args()

    config = {
        # 'input_bam': args.input_bam,
        'input_fq': args.input_fq,
        'input_fq2': args.input_fq2,
        'num_cpus': multiprocessing.cpu_count(),

        'steps': {
            'biobloomcategorizer': {
                'input_bf': args.input_bf,
                'bf_name': re.sub('\.bf$', '', os.path.basename(args.input_bf)),
            },

            'transabyss': {
                'num_reads_cutoff': args.num_reads_too_small_cutoff,
                'kmer_sizes': args.transabyss_kmer_sizes,
            },

            'align_contigs_to_genome': {
                'reference_genome_bwa_index': args.reference_genome_bwa_index
            },

        },

    }

    config['output_dir'] = os.path.join(
        os.path.dirname(config['input_fq']), '{0}-results'.format(project_id))

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
