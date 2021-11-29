#! /usr/bin/env python

import incar

inc = incar.incar().settag(
        """ialgo = 48,
        lreal = Auto,
        nelmin = 4,
        isym = 0,
        nsw = 40,
        potim = 1.0,
        tebeg = 273,
        prec = Normal,
        maxmix = 60,
        encut = 400"""
        )

print """Set
ENCUT [400]
NSW [40]
POTIM [1.0]
TEBEG [273]
MAXMIX [60]
to the values you desire."""

inc.write()

kpoints = open('KPOINTS', 'w')

kpoints.write("""Gamma
0
Gamma
1 1 1
0 0 0
"""
)
kpoints.close()
