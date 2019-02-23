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
	virtualenv -p "`which python2.7`" "$VENVDIR"
fi

source "$VENVDIR"/bin/activate
if [ "`which python2.7`" != "$VENVDIR/bin/python2.7" ]; then
	echo "ERROR: Virtualenv is not there or couldn't activate. Aborting so LAMMPS doesn't install on system python."
	echo "       Re-run this script to finish installation after sorting out whatever is wrong with virtualenv."
	exit
fi

pip install scipy numpy pyquaternion matplotlib mpi4py

echo "...done!"

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

#get the cosine/squared pair style
#git clone -b master https://github.com/Saric-Group/lammps_pair_cosine_squared.git "$LAMMPSDIR"/temp
#mv "$LAMMPSDIR"/temp/src/* "$LAMMPSDIR"/src
#rm -rf "$LAMMPSDIR"/temp

cd "$LAMMPSDIR"/src
make clean-all
make purge
make package-update
make no-all
make yes-mc yes-molecule yes-rigid yes-opt yes-misc yes-user-misc yes-python
# why does mpi not work?!?!
make -j4 serial #LMP_INC="-DLAMMPS_PNG -DLAMMPS_JPEG -DLAMMPS_FFMPEG" JPG_LIB="-lpng -ljpeg"
make -j4 serial mode=shlib LMP_INC="-DLAMMPS_EXCEPTIONS"
make install-python

echo "...done!"
deactivate
cd "$STARTDIR"
