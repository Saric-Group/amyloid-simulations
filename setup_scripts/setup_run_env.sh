#!/bin/bash -e

# This script sets up virtualenv to PROJECT_ROOT/venv and LAMMPS to
# PROJECT_ROOT/lammps-local, and installs this LAMMPS into the venv,
# with everything needed to run the programs of this project.

STARTDIR="$(realpath .)"
cd "$( dirname "${BASH_SOURCE[0]}" )/.."
BASEDIR="$(realpath .)"

# 1) setup python virtual environment

VENVDIR="$BASEDIR"/venv
if [ ! -d $VENVDIR ]; then
	echo "Installing and setting up virtualenv..."
	venv_setup=1
	python3 -m venv "$VENVDIR"
fi

source "$VENVDIR"/bin/activate
if [ "`which python3`" != "$VENVDIR/bin/python3" ]; then
	echo "`which python3` != $VENVDIR/bin/python3 !!"
	echo "ERROR: Virtualenv is not there or couldn't activate. Aborting so LAMMPS doesn't install on system python."
	echo "       Re-run this script to finish installation after sorting out whatever is wrong with virtualenv."
	cd "$STARTDIR"
	exit
fi

if [[ $venv_setup == 1 ]]; then
	pip install scipy pyquaternion mpi4py #matplotlib venv
	echo "...done!"
fi

# 2) download/update LAMMPS

LAMMPSDIR="$BASEDIR"/lammps-local
if [ ! -d $LAMMPSDIR ]; then
	echo "Fetching LAMMPS source code ..."
	git clone -b develop https://github.com/erozic/lammps.git "$LAMMPSDIR"
	cd "$LAMMPSDIR"
else
	echo "Updating LAMMPS source code ..."
	cd "$LAMMPSDIR"
	git checkout develop
	git pull origin develop
fi
echo "... done!"

# 3) build and install LAMMPS

echo "Building and installing LAMMPS ..."
BUILDDIR=""$LAMMPSDIR"/build"
mkdir "$BUILDDIR" 2>&1 || true
cd "$BUILDDIR"

FAIL="false" # for clean-up if build fails
cmake ../cmake/ || FAIL="true" 
Cmake_vars="$Cmake_vars"' -D BUILD_SHARED_LIBS=yes'
Cmake_vars="$Cmake_vars"' -D BUILD_MPI=yes'
Cmake_vars="$Cmake_vars"' -D BUILD_OMP=yes'
Cmake_vars="$Cmake_vars"' -D BUILD_TOOLS=no'
Cmake_vars="$Cmake_vars"' -D WITH_JPEG=yes'
Cmake_vars="$Cmake_vars"' -D WITH_PNG=yes'
Cmake_vars="$Cmake_vars"' -D WITH_GZIP=yes'
Cmake_vars="$Cmake_vars"' -D WITH_FFMPEG=yes'
cmake $Cmake_vars . || FAIL="true"
LAMMPS_packages='-D PKG_MOLECULE=on'
LAMMPS_packages="$LAMMPS_packages"' -D PKG_RIGID=on'
LAMMPS_packages="$LAMMPS_packages"' -D PKG_MC=on'
LAMMPS_packages="$LAMMPS_packages"' -D PKG_EXTRA-PAIR=on'
LAMMPS_packages="$LAMMPS_packages"' -D PKG_PYTHON=on'
cmake -C ../cmake/presets/all_off.cmake $LAMMPS_packages . || FAIL="true"
echo
cmake --build . --target all || FAIL="true" # problem: one core only
#make -j 4 all || FAIL="true" # problem: how many cores to use ??
echo

if [ "$FAIL" = "false" ]; then
  echo "$(tput setaf 2)$(tput smso)Build successfull!$(tput rmso) Trying to install lammps python library...$(tput sgr0)"
  echo
  make install-python || {
    echo
    echo "$(tput setaf 1)$(tput smso)ERROR:$(tput rmso) LAMMPS python lib installation failed. Check if wheel file generated in $BUILDDIR and try to install it manually.$(tput sgr0)"
  }
  echo
  echo "... done!"
else
  echo "... well, seems something went wrong during build. Check output to see what exactly."
fi

# finish
deactivate
cd "$STARTDIR"
