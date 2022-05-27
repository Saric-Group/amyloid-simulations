#!/bin/bash -e

# This script sets up virtualenv to PROJECT_ROOT/venv and LAMMPS to
# PROJECT_ROOT/lammps-local, and installs this LAMMPS into the venv,
# with everything needed to run the programs of this project.

STARTDIR="$(pwd)"
cd "$( dirname "${BASH_SOURCE[0]}" )/.."
BASEDIR="$(pwd)"

VENVDIR="$BASEDIR"/venv
if [ ! -d $VENVDIR ]; then
	echo "Installing and setting up virtualenv..."
	venv_setup=1
	virtualenv -p "`which python2.7`" "$VENVDIR"
fi

source "$VENVDIR"/bin/activate
if [ "`which python2.7`" != "$VENVDIR/bin/python2.7" ]; then
	echo "ERROR: Virtualenv is not there or couldn't activate. Aborting so LAMMPS doesn't install on system python."
	echo "       Re-run this script to finish installation after sorting out whatever is wrong with virtualenv."
	exit
fi

if [[ $venv_setup == 1 ]]; then
	pip install virtualenv scipy pyquaternion mpi4py #matplotlib
	echo "...done!"
fi

LAMMPSDIR="$BASEDIR"/lammps-local
if [ ! -d $LAMMPSDIR ]; then
	echo "Fetching, building and installing LAMMPS..."
	git clone -b develop https://github.com/erozic/lammps.git "$LAMMPSDIR"
	cd "$LAMMPSDIR"
else
	echo "Updating, rebuilding and reinstalling LAMMPS..."
	cd "$LAMMPSDIR"
	git checkout develop
	git pull origin develop
fi

cd "$LAMMPSDIR"/src
make clean-all
make purge
make package-update
make no-all
make yes-molecule yes-rigid yes-mc yes-python yes-extra-pair
#make -j4 mpi #LMP_INC="-DLAMMPS_PNG -DLAMMPS_JPEG -DLAMMPS_FFMPEG" JPG_LIB="-lpng -ljpeg"
make -j4 mpi mode=shared LMP_INC="-DLAMMPS_EXCEPTIONS"
make install-python

echo "...done!"
deactivate
cd "$STARTDIR"
