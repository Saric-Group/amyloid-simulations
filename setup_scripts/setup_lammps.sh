#!/bin/bash -e

# Will download (or just update) and (re)install the latest LAMMPS stable version
# to a desired location (in a "lammps-stable" folder).
# It will also create a "lammps-stable" link in "~/.local/bin" to the serial version
# of the installed LAMMPS.

STARTDIR="$(pwd)"

if [ -z "$1" ]; then
	echo "A location for the 'lammps-stable' folder is required as an argument!"
	exit
else
	LAMMPSDIR="$1"/lammps-stable
fi

if [ ! -d $LAMMPSDIR ]; then
	git clone -b stable https://github.com/lammps/lammps.git $LAMMPSDIR
else
	cd "$LAMMPSDIR"
	git checkout stable
	git pull
fi

cd "$LAMMPSDIR"/src
make clean-all
make purge
make package-update
make no-all
make yes-dipole yes-kspace yes-mc yes-molecule yes-rigid yes-opt yes-misc yes-user-misc yes-python
make -j4 serial #LMP_INC="-DLAMMPS_PNG -DLAMMPS_JPEG -DLAMMPS_FFMPEG" JPG_LIB="-lpng -ljpeg"
make -j4 serial mode=shlib LMP_INC="-DLAMMPS_EXCEPTIONS"
# why does mpi not work?!?!
make install-python

cd "$LAMMPSDIR"/python
python install.py #system-level python install

ln -fs "$LAMMPSDIR"/src/lmp_serial ~/.local/bin/lammps-stable

cd "$STARTDIR"
