#! /usr/bin/env python

import os, sys
from shutil import copy2 as copy
import argparse

"""
code to submit a job given the file 'ready' exists in the given
directory/directories
"""


def check_and_submit(dirs, cmd):

    # get present working directory
    pwd = os.getcwd()
    for dir in dirs:
        readyfile = os.path.join(dir, 'ready')
        if os.path.exists(readyfile):
            src = os.path.join(pwd, submitfile)
            dst = os.path.join(dir, 'sub')
            copy(src, dst)
            # change directory to dir
            os.chdir(dir)
            # do qsub submitfile
            print 'Submitting job in %s' % (dir)
            os.system(cmd + ' sub')
            # remove the 'ready' file
            os.remove('ready')
            # get back to original directory
            os.chdir(pwd)

    return None

if __name__ == '__main__':

    parser = argparse.ArgumentParser()

    parser.add_argument('-d', '--dirs', nargs='+', required=True, 
                help='directories in which to submit')
    parser.add_argument('-f', '--submitfile', default='sub', nargs='?',
                help='submit file to submit using qsub (default = "sub")')
    parser.add_argument('-c', '--command_to_use', default='qsub',
                help='command to use with the submit file (default = "qsub")')

    args = parser.parse_args()

    # directories to check for the file 'ready' taken from the command line
    dirs = args.dirs
    print dirs
    # submit file for use with qsub, to do qsub <submitfile>
    submitfile = args.submitfile
    cmd = args.command_to_use

    check_and_submit(dirs, cmd)
