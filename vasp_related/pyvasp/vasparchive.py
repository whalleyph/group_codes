import tarfile

def get_files(tgzs, file):

    """
    Returns a list of file objects of the <file> extracted from the 
    list of <tgzs>
    """

    fileobjs = []
    for tgz in tgzs:
        if tarfile.is_tarfile(tgz):
            number = arg[:-4]
            filename = file + '.' + number
            file = tarfile.open(arg).extractfile(outcarname)
            fileobjs.append(file)
        else:
            print "%s not a tar file" % (tarfile)

    return fileobjs
