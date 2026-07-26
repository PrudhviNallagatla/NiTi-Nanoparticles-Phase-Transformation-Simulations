#!/bin/bash
set -euo pipefail

export DEBIAN_FRONTEND=noninteractive
export TZ=Etc/UTC

sudo rm -rf /usr/local/bin/lmp /usr/local/bin/lammps /usr/local/lib/liblammps* /usr/local/lib64/liblammps* /tmp/lammps 2>/dev/null || true
pip3 uninstall -y lammps 2>/dev/null || true
sudo ldconfig

sudo apt-get update -qq && sudo apt-get upgrade -y -qq
sudo apt-get install -y --no-install-recommends \
    build-essential cmake git git-lfs wget curl vim nano tmux btop unzip zip pigz zstd xz-utils jq bc time \
    gfortran g++ gcc make pkg-config libopenmpi-dev openmpi-bin libblas-dev liblapack-dev libfftw3-dev \
    libpng-dev libjpeg-dev libz-dev python3 python3-dev python3-pip python3-venv libhdf5-serial-dev hdf5-tools \
    libnetcdf-dev netcdf-bin ffmpeg gnuplot-nox libkim-api-dev openkim-models

if ! command -v zellij &> /dev/null; then
    sudo snap install zellij --classic 2>/dev/null || true
fi

python3 -m pip install --break-system-packages --upgrade pip setuptools wheel > /dev/null
pip3 install --break-system-packages --no-cache-dir \
    numpy scipy matplotlib pandas ase ovito MDAnalysis pymatgen freud-analysis \
    jupyter jupyterlab tqdm pyyaml h5py netcdf4 scikit-learn lmfit atomman phonopy Cython

cd /tmp
git clone --depth 1 --branch stable_2Aug2023_update3 https://github.com/lammps/lammps.git
cd lammps
rm -rf build && mkdir -p build && cd build

export CC=mpicc CXX=mpicxx FC=mpifort
cmake ../cmake -D CMAKE_BUILD_TYPE=Release -D CMAKE_CXX_COMPILER=mpicxx -D CMAKE_CXX_FLAGS="-O3 -march=native -DNDEBUG" \
    -D BUILD_MPI=yes -D BUILD_OMP=yes -D PKG_OPENMP=yes -D PKG_INTEL=yes -D PKG_OPT=yes -D PKG_MEAM=yes -D PKG_MANYBODY=yes \
    -D PKG_PTM=yes -D PKG_VORONOI=yes -D DOWNLOAD_VORONOI=yes -D PKG_REPLICA=yes -D PKG_MC=yes -D PKG_PYTHON=yes \
    -D PKG_PLUMED=yes -D DOWNLOAD_PLUMED=yes -D PKG_IO=yes -D PKG_MPIIO=yes -D PKG_COMPRESS=yes -D PKG_EXTRA-COMPUTE=yes \
    -D PKG_EXTRA-DUMP=yes -D PKG_EXTRA-FIX=yes -D WITH_GZIP=yes -D WITH_JPEG=yes -D WITH_PNG=yes -D WITH_FFMPEG=yes

make -j"$(nproc)"
sudo make install

echo "/usr/local/lib" | sudo tee /etc/ld.so.conf.d/lammps.conf > /dev/null
sudo ldconfig

export LD_LIBRARY_PATH=/usr/local/lib:${LD_LIBRARY_PATH:-}
export LAMMPS_LIB_PATH=/usr/local/lib

cd ../python
pip3 install --break-system-packages . || { cd ../build; sudo make install-python; }

if ! grep -q "LAMMPS_ULTIMATE" ~/.bashrc; then
cat >> ~/.bashrc << 'EOF'
export LD_LIBRARY_PATH=/usr/local/lib:$LD_LIBRARY_PATH
export PYTHONPATH=$PYTHONPATH:/usr/local/lib/python3/dist-packages
alias lmp='lmp'
alias lmp_mpi='mpirun -np $(nproc) lmp'
alias lmp_intel='mpirun -np $(nproc) lmp -sf intel -pk intel 0 mode mixed'
EOF
fi

cd ~
sudo rm -rf /tmp/lammps