#!/bin/bash -e

# Will install "virtualenv" to a desired location (in a "venv" folder).
# The virtualenv will have scipy, numpy, matplotlib and mpi4py.
# The script assumes there is a "lammps-stable" exe link on the system which
# points a LAMMPS (serial) application in a "src" folder of a LAMMPS installation.

STARTDIR="$(pwd)"

if [ -z "$1" ]; then
	echo "A location for the 'venv' folder is required as an argument!"
	exit
fi

cd "$1"
VENVDIR="`pwd`"/venv

if [ ! -d $VENVDIR ]; then
	virtualenv -p "`which python2.7`" "$VENVDIR"
fi

source "$VENVDIR"/bin/activate
if [ "`which python2.7`" != "$VENVDIR/bin/python2.7" ]; then
	echo "ERROR: Virtualenv is not there or couldn't activate. Aborting so LAMMPS doesn't install on system python. Re-run this script to finish installation after sorting out whatever is wrong with virtualenv."
	exit
fi

pip install scipy numpy pyquaternion matplotlib mpi4py

LAMMPSDIR="$(dirname "$(readlink -f "`which lammps-stable`")")"/..

cd "$LAMMPSDIR"/python
python install.py

deactivate
cd "$STARTDIR"
