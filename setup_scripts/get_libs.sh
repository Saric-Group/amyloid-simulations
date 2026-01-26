#!/bin/bash -e

# Will download the latest necessary libraries to the "libs" folder in
# the project root folder. These are currently:
#  - lammps_multistate_rods

STARTDIR="$(pwd)"
cd "$( dirname "${BASH_SOURCE[0]}" )/.."
BASEDIR="$(pwd)"

# 1) get "lammps multistate rods" library (put it in virtual env)

py_version=$(python --version)
[[ $py_version =~ Python\ ([[:digit:]]+)\.([[:digit:]]+)\.([[:digit:]]+) ]] && py_version="${BASH_REMATCH[1]}.${BASH_REMATCH[2]}"

LIBDIR="$BASEDIR/venv/lib/python${py_version}/site-packages"
mkdir -p "$LIBDIR"

git clone -b develop https://github.com/Saric-Group/lammps_multistate_rods.git "$LIBDIR/temp"

rm -rf "$LIBDIR/lammps_multistate_rods"
mv "$LIBDIR/temp/lammps_multistate_rods" "$LIBDIR/lammps_multistate_rods"
rm -rf "$LIBDIR/temp"

# 2) ...

# finish
cd "$STARTDIR"
