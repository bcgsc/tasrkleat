import os
import datetime
import gzip
import subprocess
import time
import select
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
        logger.info("Time spent on {0}: {1:.2f}s".format(f.__name__, et - bt))
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


def execute(cmd, flag_file=None, msg_id='', debug=False):
    """
    # http://stackoverflow.com/questions/1606795/catching-stdout-in-realtime-from-subprocess
    :param cmd: should never inlcude pipe or redirection, which would requires
    a new shell process

    This execute logs all stdout and stderr, which could look funny, especially
    when it comes to tools like aspc and wget
    """
    logger.info('executing: {0}'.format(cmd))
    # todo: should check whether cmdsp includes pipe or redirection here

    if debug:
        return
    try:
        proc = subprocess.Popen(
            cmd,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            shell=True,
            executable='/bin/bash',
            # https://docs.python.org/3.5/library/subprocess.html#subprocess.Popen.communicate
            universal_newlines=True
        )
        ioselect(proc)
        returncode = proc.returncode
        msg = '{0}: returncode: {1}. CMD: "{2}"'.format(msg_id, returncode, cmd)

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


def ioselect(proc):
    """
    select in the context io completion,
    https://docs.python.org/2/library/select.html
    """
    while True:
        ret = select.select([proc.stdout.fileno(),
                             proc.stderr.fileno()],
                            [], [])
        # print(ret)
        for fd in ret[0]:
            if fd == proc.stdout.fileno():
                line = proc.stdout.readline()
                if line:
                    # rstrip: remove newline character since logger will
                    # add one automatically
                    logger.info('stdout: ' + line.rstrip())
            if fd == proc.stderr.fileno():
                line = proc.stderr.readline()
                if line:
                    logger.warning('stderr: ' + line.rstrip())

        # check if child process has terminated
        # https://docs.python.org/3/library/subprocess.html#subprocess.Popen.poll
        if proc.poll() != None:
            # flush everything first, otherwise, some stdout/stderr may not
            # be able to be written to screen (e.g. cutadapt)
            for fd in ret[0]:
                if fd == proc.stdout.fileno():
                    for line in proc.stdout.readlines():
                        if line:
                            # rstrip: remove newline character since logger will
                            # add one automatically
                            logger.info('stdout: ' + line.rstrip())
                if fd == proc.stderr.fileno():
                    for line in proc.stderr.readlines():
                        if line:
                            logger.info('stderr: ' + line.rstrip())
            break


def fastq_too_small(fq, num_reads=50):
    num = num_reads * 4
    count = 0
    with open(fq) as inf:
        for line in inf:
            count += 1
            if count >= num:
                return False, None
    return True, count / 4
