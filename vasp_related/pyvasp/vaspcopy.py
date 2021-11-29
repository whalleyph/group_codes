#! /usr/bin/env python

import os
from shutil import copy

def copy_start(src_dir, dst_dir):
    # Copy INCAR KPOINTS POSCAR and POTCAR
    files = 'INCAR KPOINTS POSCAR POTCAR'.split()
    for file in files:
        src = os.path.join(src_dir, file)
        dst = os.path.join(dst_dir, file)
        print 'Copying %s to %s' % (src, dst)
        copy(src, dst)

    return None


def copy_final(src_dir, dst_dir, link=False):
    def copy_final_1():
        # Copy CONTCAR to POSCAR if CONTCAR exists; else copy POSCAR
        src = os.path.join(src_dir, 'CONTCAR')
        if not os.path.exists(src):
            src = os.path.join(src_dir, 'POSCAR')
        dst = os.path.join(dst_dir, 'POSCAR')
        print 'Copying %s to %s' % (src, dst)
        copy(src, dst)
        # Copy INCAR KPOINTS and POTCAR
        files = 'INCAR KPOINTS POTCAR'.split()
        for file in files:
            src = os.path.join(src_dir, file)
            dst = os.path.join(dst_dir, file)
            print 'Copying %s to %s' % (src, dst)
            copy(src, dst)
    def create_links_wavecar_chgcar():
        files = 'CHGCAR WAVECAR'.split()
        os.chdir(dst_dir)
        for file in files:
            src = os.path.join('../bulk', file)
            lnk = file
            print 'Creating symliks to %s' % (lnk)
            os.symlink(src, lnk)

    if link:
        copy_final_1()
        create_links_wavecar_chgcar()
    else:
        copy_final_1()

    return None


if __name__ == '__Main__':
    pass
