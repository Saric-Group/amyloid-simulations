#!/bin/bash -e

# Will download the latest necessary libs to the "libs" folder in the project root folder.

STARTDIR="$(pwd)"
cd "$( dirname "${BASH_SOURCE[0]}" )/.."
BASEDIR="$(pwd)"

mkdir -p "libs"
LIBSDIR="$BASEDIR/libs"

git clone -b master https://github.com/Saric-Group/lammps_multistate_rods.git "$LIBSDIR/temp"

rm -rf "$LIBSDIR/lammps_multistate_rods"
mv "$LIBSDIR/temp/lammps_multistate_rods" "$LIBSDIR/lammps_multistate_rods"
rm -rf "$LIBSDIR/temp"

cd "$STARTDIR"
