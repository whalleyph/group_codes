#! /usr/bin/env python

import eigenval

f = '/data/chaitanya/work/graphene/band/EIGENVAL'

e = eigenval.eigenval(f)

print e.parse()
