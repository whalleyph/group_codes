# Copy a file or list of files to all directories one level deep from the
# the required library. So,
# recursivecopy(filelist, dir) will copy filelist to all directories
# dir/{dir1,dir2,dir3,...}

import sys
import os
from shutil import copy2

def recursivecopy(filelist, targetdir):
    dirs = [i for i in os.listdir(targetdir) if os.path.isdir(i)]
    for dir in dirs:
        dstdir = os.path.join(targetdir, dir)
        for file in filelist:
            src = file
            dst = os.path.join(dstdir, file)
            copy2(src, dst)

    return None

if __name__ == '__main__':
    targetdir = sys.argv[1]
    filelist = sys.argv[2].split()
    recursivecopy(filelist, targetdir)
