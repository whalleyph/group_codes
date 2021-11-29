#! /usr/bin/env python

import outcar as out
import argparse

if __name__ == '__main__':
    parser = argparse.ArgumentParser()

    parser.add_argument('-i', '--input_file', default='OUTCAR')

    args = parser.parse_args()
    outcar = out.outcar().read(args.input_file)

    print outcar.get_forces()
