#!/bin/bash -e

# Will download (or just update) and (re)install the latest private LAMMPS fork
# to a desired location (in a "lammps-develop" folder).
# It will build the serial version, both the executable and the library.
# It will also create a "lammps-develop" link in "~/.local/bin" to the serial version
# of the installed LAMMPS.

STARTDIR="`pwd`"

if [ -z "$1" ]; then
	echo "A location for the 'lammps-develop' folder is required as an argument!"
	exit
fi

cd "$1"
LAMMPSDIR="`pwd`"/lammps-develop

if [ ! -d $LAMMPSDIR ]; then
	git clone -b develop https://github.com/erozic/lammps.git "$LAMMPSDIR"
	cd "$LAMMPSDIR"
else
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

cd "$LAMMPSDIR"/python
python install.py #system-level python install

ln -fs "$LAMMPSDIR"/src/lmp_serial ~/.local/bin/lammps-develop

cd "$STARTDIR"
