import sys
import os
import subprocess as sp
import time
import paramiko
from scp import SCPClient
import logging
import argparse as ap
from time import sleep, time
import random
import string


def submit(ssh, scp, files_to_send, wd, x):
    files_to_send.append(x)
    _, pwd, _ = ssh.exec_command("pwd")
    _ = ssh.exec_command(' '.join(['mkdir -p', wd.rstrip("/")]))
    scp.put(files_to_send, remote_path=wd)
    i, o, e = ssh.exec_command(' '.join(['cd', wd, '&&', 'sbatch', x]))
    _ = ssh.exec_command("cd "+pwd.readlines()[0].strip())
    out = o.readlines()
    id = out[-1].strip().split()[-1]
    if id.isnumeric():
        return id
    else:
        print "SBATCH filed. Exitting..."
        exit(1)


def check(ssh, jid, period, r):
    while True:
        utime = time()
        current = "running"+r
        os.system('echo '+str(utime)+" > "+current)
        i, o, e = ssh.exec_command("squeue -u $USER")
        out = o.readlines()
        ids = []
        for line in out:
            l = line.strip().split()
            ids.append(l[0])
        if jid not in ids:
            break
        sleep(period)
    os.system("rm "+current)
    return True


def collect(scp, wd, files_to_get):
    pwd = wd.rstrip("/")
    for f in files_to_get:
        fs = pwd+"/"+f
        scp.get(fs, local_path='.')
    pass

if __name__ == "__main__":
    parser = ap.ArgumentParser(
        description="Submit jobs and recieve the results with python! The script will remain running unless: 1. The job stops serverside; 2. The server conncection is broken; 3. You kill it."
    )
    parser.add_argument(
        "-s",
        "--server",
        help="Server address",
        required=True,
    )

    parser.add_argument(
        "-k",
        "--key",
        help="The absolute location of your private key giving access to the server",
        required=True,
        type=str
    )

    parser.add_argument(
        "-u",
        "--user",
        help="",
        required=True,
    )

    parser.add_argument(
        "-d",
        "--directory",
        help="The absolute path of your working directory on the server",
        required=True,
        type=str,
    )

    parser.add_argument(
        "-x",
        "--submission",
        help="The submission file for the task scheduler",
        required=True,
    )

    parser.add_argument(
        "-f",
        "--files",
        help="Files to be copied over",
        required=False,
        nargs="+",
    )

    parser.add_argument(
        "-t",
        "--period",
        help="The time period (in seconds) in which to repeatedly check the server for completion. Def:. 30",
        default=30,
        type=int,
        required=False,
    )

    parser.add_argument(
        "-o",
        "--output",
        help="List of output files to be collected from the server once the job is complete.",
        required=True,
        nargs="+",

    )  # We need str to avoid picking up a local file list from use of wildcards

    args = parser.parse_args()
    # print unicode("/home/chem/msrzvr/.ssh/avon_rsa2")
    pk = unicode(args.key.strip())

    ssh = paramiko.SSHClient()
    k = paramiko.RSAKey.from_private_key_file(
    pk,
    )

    pool = string.ascii_lowercase+string.ascii_uppercase
    rnd = ''.join(random.choice(pool) for i in range(8))

    ssh.set_missing_host_key_policy(paramiko.AutoAddPolicy())
    ssh.connect(hostname=args.server, username=args.user, pkey=k)
    scp = SCPClient(ssh.get_transport())
    dir = args.directory.rstrip("/")+"/calc_"+rnd+"/"
    jid = submit(ssh, scp, files_to_send=args.files, wd=dir, x=args.submission,)
    done = check(ssh, jid, args.period, rnd)
    if done:
        collect(scp, dir, args.output)

"""
Calypso needs OUTCAR and CONTCAR
"""
