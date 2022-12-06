"""
parallel_mcsqs.py by Ben Saunders 2022

DESCRIPTION
    Uses concurrent futures to call n parallel processes of ATAT's mcsqs program. In
    principle, one only needs to change the comment in the function 'mcsqs' to run
    any command in parallel.

USAGE
    python parallel_mcsqs.py -n <number of sites required by mcsqs> -p <number of parallel processes>
    
    (In future, script could read 'reminder.info' to obtain minimum number of sites to avoid mcsqs
    throwing errors?)

"""
import argparse as ap
import subprocess as sp
import concurrent.futures as cf
from time import sleep


def run(command):
    sleep(0.1)  # Just in case
    p = sp.Popen(command, stdout=sp.PIPE, encoding="utf8")
    output, err = p.communicate()
    return output, err


def mcsqs(n, id):
    cmd = ["mcsqs", f"-n={n}", f"-ip={id}"]
    print(f"MCSQS PARALLEL - Process ID {id}\t... Running command: {cmd}")
    o, e = run(command=cmd)
    print(f"MCSQS PARALLEL - Process ID {id} has been terminated.")
    return (o, e)


def paralellise(sites, ids, workers=1):
    arg_list = [[sites, id] for id in ids]

    with cf.ProcessPoolExecutor(max_workers=workers) as executor:
        x = executor.map(mcsqs, *zip(*arg_list))
    return list(x)


if __name__ == "__main__":
    parser = ap.ArgumentParser(description="Set multiple mcsqs sessions off at once.")

    parser.add_argument(
        "-p",
        "--procs",
        type=int,
        help="Number of sessions of mcsqs to start",
        default=1,
    )
    parser.add_argument(
        "-n", "--sites", type=int, help="Size of the bestsqs cell.", default=0
    )
    args = parser.parse_args()

    process_ids = [*range(0, args.procs, 1)]
    paralellise(args.sites, process_ids, args.procs)
