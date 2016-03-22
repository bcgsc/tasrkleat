import os
import datetime
import gzip
import subprocess
import time
from functools import update_wrapper

import logging
logger = logging.getLogger(__name__)


def decorator(d):
    "Make function d a decorator: d wraps a function fn."
    def _d(fn):
        return update_wrapper(d(fn), fn)
    update_wrapper(_d, d)
    return _d


@decorator
def timeit(f):
    """time a function, used as decorator"""
    def new_f(*args, **kwargs):
        bt = time.time()
        r = f(*args, **kwargs)
        et = time.time()
        delta_t = datetime.timedelta(seconds=(et - bt))
        logger.info("Time spent on {0}: {1}".format(f.__name__, delta_t))
        return r
    return new_f


def gzip_compress(input_file):
    """name as such to differentiate from gzip module"""
    with gzip.open(input_file + '.gz', 'wb') as opf:
        with open(input_file, 'rb') as inf:
            opf.write(inf.read())
    os.remove(input_file)


def touch(fname, cmd=None, times=None):
    """
    Similar to nix command, touch, to create an empty file, but also added some
    meta data to the touched file
    """
    with open(fname, 'a') as opf:
        opf.write('created: {0}\n'.format(datetime.datetime.now()))
        opf.write('location of code execution: {0}\n'.format(
            os.path.abspath('.')))
        if cmd:
            opf.write('{0}\n'.format(cmd))
        os.utime(fname, times)


def execute(cmd, flag_file=None, msg_id='#', debug=False):
    """
    This execute doesn't log all stdout, which could look funny, especially
    when it comes to tools like aspc and wget
    """
    logger.info('{0}: executing: {1}'.format(msg_id, cmd))
    if debug:                   # only print out cmd
        return
    try:
        returncode = subprocess.call(cmd, shell=True, executable="/bin/bash")
        msg = '{0}: returncode: {1}; CMD: "{2}"'.format(msg_id, returncode, cmd)
        if returncode != 0:
            logger.error(msg)
        else:
            logger.info(msg)
            if flag_file is not None:
                touch(flag_file, cmd)
        return returncode
    except OSError as err:
        logger.exception(
            '{0}: failed to start, raising OSError {1}. '
            'CMD: "{2}"'.format(msg_id, err, cmd))
