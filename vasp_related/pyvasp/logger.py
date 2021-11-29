#! /usr/bin/env python

from datetime import datetime

class Logger:
    def __init__(self, logfile):
        self.log = logfile
        open(self.log, 'w').write('Starting\n')

    def write(self, string):
        string = datetime.today().isoformat().split('.')[0] + ': ' + string
        open(self.log, 'a').write(string.strip()+'\n')
