# encoding: utf-8
'''
Takes a LAMMPS dump file and keeps only the good part of it (from start to where something went wrong)

Created on 31 March 2019

@author: Eugen Rožić
'''

import os
import argparse
from lammps_multistate_rods.tools import parse_dump_file, write_dump_snapshot

parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument('in_files', nargs='+',
                    help='path(s) of the dump file(s) to clean')
args = parser.parse_args()

for dump_filename in args.in_files():
    try:
        for parse_out in parse_dump_file(dump_filename):
            write_dump_snapshot(parse_out, dump_filename + '_', append = True)
    except Exception as e:
        print e
        print "{:s} was good until timestep {:d}".format(dump_filename, parse_out[0])

    os.remove(dump_filename)
    os.rename(dump_filename + '_', dump_filename)
