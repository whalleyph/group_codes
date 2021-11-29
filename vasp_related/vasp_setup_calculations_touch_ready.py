#!/usr/bin/env python
# encoding: utf-8

import os

for i in os.walk('.'):
    if not i[1]:
        open(os.path.join(i[0], 'ready'), 'w')
