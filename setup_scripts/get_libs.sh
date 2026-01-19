#!/bin/bash -e

# Will download the latest necessary libraries to the "libs" folder in
# the project root folder. These are currently:
#  - lammps_multistate_rods

STARTDIR="$(pwd)"
cd "$( dirname "${BASH_SOURCE[0]}" )/.."
BASEDIR="$(pwd)"

# 1) get "lammps multistate rods" library (put it in virtual env)

LIBDIR="$BASEDIR/venv/lib/python3.12/site-packages"
mkdir -p "$LIBDIR"
rm -rf "$LIBDIR/lammps_multistate_rods"

git clone -b develop https://github.com/Saric-Group/lammps_multistate_rods.git "$LIBDIR/lammps_multistate_rods"

# 2) ...

# finish
cd "$STARTDIR"
