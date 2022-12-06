"""
=============================
=|| CALYPSO REMOTE RUNNER ||=
=============================

calypso_remote_runner.py
Benedict Saunders, 2022

Script to perform CALYSPO in split mode, running all QM calcultions on a remote cluster over SSH,
using an ssh key file (via the paramiko module).

How it works:
To begin, CALYPSO needs to be in split mode. This is checked for in the script, and will terminate if it is not on.
You also need to provide the POTCAR (ordered the same as the elements in your input.dat) and the INCAR files at
various accuracies, indexed according to their accuracy.

After checking for the inputs, CALYPSO is called which will generated a number of POSCAR files. Each of these is copied
to the remote server, along with the POTCAR, INCARs and submission file, into a new directory, labelled with the run ID,
generation number and structure number. Each of these vasp runs are queued using the supplied submission file, and the
SLURM job ID is recorded.
The command squeue is run periodically on the remote server, and the list of job IDs queued or running is checked to see
if they match the job IDs of the above jobs. If any do match, the script continues checking; if not, the jobs are deemed
complete, and CONTCAR and OUTCAR are downloaded from each directory, and renamed according to the structure number.
Once the indexed OUTCARs and CONTCARs have been collected, CALYPSO is called again and the process repeats until CALYPSO
deems the output to be the energy minimum.

"""

VERSION = "0.1.9"

import os
import subprocess as sp
import paramiko
from scp import SCPClient
import argparse as ap
from time import sleep, strftime, time
import random
import string
import glob
from pyfiglet import Figlet

def get_calypso_settings():
    print "Reading CALYPSO inputs...\n"
    split_present = False
    with open('input.dat', 'r') as input_dat:
        lines = input_dat.readlines()
    for line in lines:
        l = line.strip()
        if "plit" in l:
            split_mode= l.split("=")[-1].strip()
            split_present = True
            if split_mode != "T":
                print "Enable split mode in your input.dat"
            break
    if not split_present:
        print "Could not find the split mode switch in your input.dat"
    try:
        with open("step",'r') as step:
            curr = step.readlines()
        current = curr[0].strip()
        print "Continuing as generation "+current
    except:
        current = 1
        print "CALYPSO step file not found. Assuming generation 1"
    return {"gen":int(current),}

def get_max_accuracy():
    incars = glob.glob("INCAR_*")
    if len(incars) < 1:
	    print "No INCAR files found. Quitting..."
	    exit(1)
    accs = [int(inc.replace("INCAR_", "")) for inc in incars]
    return max(accs), incars

def print_gen(msg, gen, indent=1):
    print "G"+str(gen)+" >>"+" "*indent+str(msg)

def submit(ssh, scp, files_to_send, wd, x):
    """
    This function creates a directory (or nested, hence the -p) before using SCP to put the input files
    there, and submitting the job in said directory. There are (hopefully now 'were') a few issues on
    the remote side of things, hence the sleeps and remote 'ls' call.
    """
    files_to_send.append(x)
    sleep(0.1)
    mi, mo, me = ssh.exec_command(' '.join(['mkdir -vp', wd.rstrip("/")]))
    mout = mo.readlines()
    for m in mout:
        if not "mkdir: created directory " in m.strip():
            raise Exception("FATAL: Remote directory not created.")

    ################################################################################################################
    li, lo, le = ssh.exec_command('ls calypso_remote')
    """
    There seems to be a bug which occurs between the remote mkdir and scp.put which only seems to be alliviated (?) by
    calling another command on the server in between. I went for listing the parent directory as we can then
    retrieve the list if we need to check it out at some point
    """
    ################################################################################################################
    sleep(1)
    scp.put(files_to_send, remote_path=wd)
    sleep(1)
    i, o, e = ssh.exec_command(' '.join(['cd', wd, '&&', 'sbatch', x]))
    out = o.readlines()
    id = out[-1].strip().split()[-1]
    if id.isnumeric():
        return id
    else:
        print "SBATCH failed. Exitting..."
        exit(1)


def check(ssh, submitted_jids, sgen, t, dirs):
    """
    Uses SQUEUE over SSH to retrieve the Job IDs of any running jobs, and compares them to the list
    obtained from when the jobs were submitted (or jobs.current if check is called from
    collect_previous()).
    Once the retrieved job IDs no longer appear in the SQUEUE list, the function breaks out of the
    dangerous while:True loop.
    """
    nums = []
    for dir in dirs:
        dirsplit = dir.rstrip("/").split("_")[-1]
        nums.append(dirsplit.replace("struct", "").strip())
    reported = []
    completed = []
    while True:
        i, o, e = ssh.exec_command("squeue -u $USER")
        timestamp = str(strftime('%H:%M:%S'))
        out = o.readlines()

        # Writing output of squeue to file
        ids = []
        with open('squeue.last', 'w') as sqlast:
            sqlast.write('SQUEUE at '+timestamp+"\n\n")
            for line in out:
                sqlast.write(line)
                l = line.strip().split()
                ids.append(l[0])
        
        # Iterating though the first column of the SQUEUE output
        for jid in submitted_jids:
            if jid not in ids and jid not in completed:
                completed.append(jid)
        if len(completed) <= len(submitted_jids):
            for jid in completed:
                if jid not in reported:
                    reported.append(jid)
                    sid = nums[submitted_jids.index(jid)]
                    print_gen("Structure "+str(sid)+" completed ("+str(reported.index(jid)+1)+"/"+str(len(submitted_jids))+")", sgen, indent=3)
        if len(completed) == len(submitted_jids):
            break
        sleep(t)
    return completed

def collect(scp, wds, files_to_get, sgen):

    max_acc = get_max_accuracy()
    failures_list = []
    """
    As the name suggests, here we collect the files in files_to_collect from the parent directory
    of the calculation, and appent _NUMBER to the filename, depending on the number CALYPSO gave
    the initally POSCAR.
    """
    ### DIRNAME FORMAT -> defined_dir/calypso_run_<ID>_gen<GEN>_struct<NUM>/"
    ### We can only the below concisely because the files_to_get locations are defined in the bash submission

    print_gen("Jobs finished. Getting outputs: "+', '.join(files_to_get), sgen)
    for dir in wds:
        dirsplit = dir.rstrip("/").split("_")[-1]
        num = dirsplit.replace("struct", "").strip()
	print_gen("Structure "+str(num)+":", sgen, indent=1)
        pwd = dir.rstrip("/")
        for f in files_to_get:
            fs = pwd+"/"+f
            fls = fs.lstrip()
            try:
                scp.get(fls)
                print_gen("Retrieved "+fls, sgen, indent=3)
            except:
                print_gen("Retrieval of "+fls+" failed.", sgen, indent=3)
                print_gen("Attempting directly from the calculation directory (Job might have been killed).", sgen, indent=3)  
                from_calc_dir = pwd+"/acc_"+str(max_acc)+"/"+f
                fls_c = from_calc_dir.lstrip()
                try:
                    scp.get(fls_c)
                except:
                    print_gen("Direct retrieval of "+fls+" failed. Skipping...", sgen, indent=5)
                    failures_list.append(int(num))
            sleep(0.1)
            os.system('mv '+f+" "+f+"_"+str(int(num)))
    failures = list(set(failures_list))
    fails = len(failures)
    print_gen("Files downloaded. "+str(fails)+" failures:", sgen)
    if len(failures) > 0:
        failures.sort()
        for idx, s in enumerate(failures):
            print_gen("Structure "+str(s).rjust(2, "0"), sgen) 


def collect_previous(ssh, scp, sgen, t, files_to_get):
    """
    For whatever reason, the main loop of the script may stop, e.g. computer restart or accidental killing
    of the job etc. Using the file 'jobs.current', we can see the location of various jobs and their job
    IDs from before the script was previously killed. 
    Like the check function above, we do an initial real of SQUEUE for any matching job IDs. If some
    previous jobs are still runnning, 'check()' is called.
    To make sure, we still have the initial POSCAR files, files_to_get is appended accordingly and
    collect is called. collect_previous() returns the run ID string to maintain consistency, and it will
    be assigned as the runID should a run need manual continuation.
    """
    prev_tmp = int(sgen)
    prev = str(prev_tmp-1)

    print_gen("Checking for currently running jobs from previous CALYPSO run.", sgen)
    jids, dirs = [], []
    with open("jobs.current", 'r') as jc_r:
        lines = jc_r.readlines()
    for l in lines:
        jid, dir = l.strip().split(":")
        jids.append(jid)
        dirs.append(dir)
    sleep(5)
    i, o, e = ssh.exec_command("squeue -u $USER")
    out = o.readlines()
    ids = []
    for sqline in out:
        ids.append(sqline.split()[0])
    if any(i_d in jids for i_d in ids):
        print_gen("Jobs still running from generation "+prev+". Checking job statuses.", sgen)
        completed = check(ssh, jids, prev, t, dirs)
    else:
        print_gen("Jobs complete from previous run.", prev)
    if "POSCAR" not in files_to_get:
        files_to_get.append("POSCAR")

    # Getting runID for continuation.
    split = dirs[-1].split("_")
    run_marker = split.index("run")+1
    prev_runID = split[run_marker]

    collect(scp, dirs, files_to_get, prev)
    return prev_runID


### BEGIN MAIN ###

if __name__ == "__main__":

    fglt = Figlet(font="larry3d")
    print fglt.renderText("  CALYPSO")
    print fglt.renderText("     REMOTE")
    print fglt.renderText("        RUNNER")
    print "CALYPSO Remote Runner, Version "+VERSION
    print "Benedict Saunders, July 2022\n"
    max_acc, globbed_incars = get_max_accuracy()
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
        "-pwd",
        "--pwd",
        help="The passphrase for your shhkey file. Who doesn't like typing out passwords in plaintext?",
        required=False,
        default=None,
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
        required=False,
        default="x.bash"
    )

    parser.add_argument(
        "-f",
        "--files",
        help="Files to be copied over",
        required=False,
        nargs="+",
        default=["POSCAR"]+globbed_incars
    )

    parser.add_argument(
        "-t",
        "--period",
        help="The time period (in seconds) in which to repeatedly check the server for completion. Def:. 30",
        default=300,
        type=int,
        required=False,
    )

    parser.add_argument(
        "-o",
        "--output",
        help="List of output files to be collected from the server once the job is complete.",
        required=False,
        nargs="+",
        default=["OUTCAR", "CONTCAR"]

    )  # We need str to avoid picking up a local file list from use of wildcards
    parser.add_argument(
        "-g",
        "--generation",
        help="Start the generation counter at something other than 1.",
        required=False,
        default=None,
    )
    parser.add_argument(
        "-continue",
        "--continue",
        dest="cont",
        help="Continue the run if this python script is killed for whatever reason. It will check and fetch strcutures from the remote server, before proceeding with a CALYPSO call.",
        required=False,
        action="store_true"
    )

    args = parser.parse_args()

    calypo_settings = get_calypso_settings()
    print "Found "+str(max_acc)+" INCAR files"

    pk = unicode(args.key.strip())
    print '\nSetting up SSH...'
    ssh = paramiko.SSHClient()
    k = paramiko.RSAKey.from_private_key_file(
        pk,
        password=args.pwd
    )
    
    ssh.set_missing_host_key_policy(paramiko.AutoAddPolicy())
    ssh.connect(hostname=args.server, username=args.user, pkey=k)
    scp = SCPClient(ssh.get_transport())
    print 'SSH connected.\n'
    
    if args.generation is not None:
        gen = int(args.generation)
    else:
        gen = int(calypo_settings["gen"])

    pool = string.ascii_lowercase+string.ascii_uppercase
    if args.cont and gen > 1:
        runID = collect_previous(ssh, scp, str(gen), args.period, args.output)
    else:
        runID = ''.join(random.choice(pool) for i in range(12))

    print 'Entering main loop at generation '+str(gen)+'...'
    while True:
        sgen = str(gen)
        sleep(5)
        print ""
        print "+ ==== GENERATION "+str(sgen)+" ==== +"
        print ""
        print_gen("Calling CALYPSO", sgen)
        sp.call('calypso.x')
        sleep(1)
        poscars = glob.glob('POSCAR_*')
        if len(poscars) < 1:
            break
        jids = []
        dirs = []
        nums = []
        sposcars = sorted(poscars)
        digits = len(list(str(len(sposcars))))
        for poscar in sposcars:
            num = int(poscar.replace("POSCAR_", ""))
            print_gen("Structure: "+str(num).rjust(digits, "0"), sgen, indent=3)
            dir = args.directory.rstrip("/")+"/calypso_run_"+runID+"_gen"+sgen.rjust(3, "0")+"_struct"+str(num).rjust(digits, "0")+"/"
            os.system('cp '+poscar+'  POSCAR')
            sleep(1)
            print_gen("Job working directory: "+dir, sgen, indent=5)
            jid = submit(ssh, scp, files_to_send=args.files, wd=dir, x=args.submission,)
            print_gen("Submitted as JobID "+str(jid)+"", sgen, indent=5)
            nums.append(num)
            dirs.append(dir)
            jids.append(jid)
        with open('jobs.current', 'w') as jc_w:
            for jid, dir in zip(jids, dirs):
                jc_w.write(str(jid)+":\t"+dir+"\n")
        os.system("rm POSCAR")
        print_gen("All submitted. Will check 'squeue' every "+str(args.period)+" second(s)", sgen)
        print_gen("Remote calculation information written to jobs.current", sgen)

        ### CHECK ###
        done = check(ssh, jids, sgen, args.period, dirs)
        print_gen("All jobs from generation completed. Collecting from server...", sgen)
        collect(scp, dirs, args.output, sgen)
	gen += 1
    print 'BROKEN OUT OF MAIN LOOP. RUN FINISHED AFTER '+sgen+" GENERATION(S)."

"""
Calypso needs OUTCAR and CONTCAR
"""
