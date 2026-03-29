#!/bin/bash
# filepath: /home/rimuru/workspace/install_lammps.sh

set -euo pipefail

# --- COLORS ---
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
BLUE='\033[0;34m'
NC='\033[0m'

log() { echo -e "${GREEN}[$(date +'%H:%M:%S')]${NC} $1"; }
warn() { echo -e "${YELLOW}[WARNING]${NC} $1"; }
error() { echo -e "${RED}[ERROR]${NC} $1"; exit 1; }

NPROC=$(nproc)
log "Detected ${NPROC} CPU cores. We will use them all."

# Set non-interactive
export DEBIAN_FRONTEND=noninteractive
export TZ=Etc/UTC

#==============================================================================
# STEP 1: PURGE OLD INSTALLATION
#==============================================================================
log "[1/7] Purging old LAMMPS installation (Scorched Earth Policy)..."

# Nuke binaries
sudo rm -f /usr/local/bin/lmp /usr/local/bin/lammps 2>/dev/null || true
# Nuke libraries
sudo rm -rf /usr/local/lib/liblammps* 2>/dev/null || true
sudo rm -rf /usr/local/lib64/liblammps* 2>/dev/null || true
# Nuke python bindings
pip3 uninstall -y lammps 2>/dev/null || true
# Nuke build dir
sudo rm -rf /tmp/lammps 2>/dev/null || true
# Clean ldconfig cache
sudo ldconfig

log "✓ Previous installation vaporized."

#==============================================================================
# STEP 2: DEPENDENCIES (UNIFIED)
#==============================================================================
log "[2/7] Updating system and installing build tools..."

sudo apt-get update -qq
sudo apt-get upgrade -y -qq

# Combined dependency list. No duplicates.
sudo apt-get install -y --no-install-recommends \
    build-essential cmake git git-lfs wget curl vim nano tmux \
    btop unzip zip pigz zstd xz-utils jq bc time \
    gfortran g++ gcc make pkg-config \
    libopenmpi-dev openmpi-bin \
    libblas-dev liblapack-dev libfftw3-dev \
    libpng-dev libjpeg-dev libz-dev \
    python3 python3-dev python3-pip python3-venv \
    libhdf5-serial-dev hdf5-tools \
    libnetcdf-dev netcdf-bin \
    ffmpeg gnuplot-nox \
    libkim-api-dev openkim-models

# Optional: Zellij
if ! command -v zellij &> /dev/null; then
    sudo snap install zellij --classic 2>/dev/null || warn "Snap not available, skipping Zellij."
fi

log "✓ Dependencies acquired."

#==============================================================================
# STEP 3: PYTHON SCIENTIFIC STACK
#==============================================================================
log "[3/7] Installing Python scientific stack..."

# Upgrade pip first
python3 -m pip install --break-system-packages --upgrade pip setuptools wheel > /dev/null

# Install the stack
# Note: Added 'wheel' and 'cython' which are often implicit deps for these packages
pip3 install --break-system-packages --no-cache-dir \
    numpy scipy matplotlib pandas \
    ase ovito MDAnalysis pymatgen freud-analysis \
    jupyter jupyterlab tqdm pyyaml h5py netcdf4 \
    scikit-learn lmfit atomman phonopy Cython

log "✓ Python stack ready."

log "[5/8] Atomsk installation SKIPPED - requires manual installation"
warn "═══════════════════════════════════════════════════════════════"
warn "  ⚠️  ATOMSK MUST BE INSTALLED MANUALLY FOR POLYCRYSTALS  ⚠️"
warn "═══════════════════════════════════════════════════════════════"
warn ""
warn "Quick install (run after this script completes):"
warn ""
warn "   git clone https://github.com/pierrehirel/atomsk.git"
warn "   cd atomsk/src"
warn "   make atomsk"
warn "   sudo make install"
warn ""
warn "Official docs: https://atomsk.univ-lille.fr/doc/en/"
warn "═══════════════════════════════════════════════════════════════"
warn ""

#==============================================================================
# STEP 4: DOWNLOAD SOURCE
#==============================================================================
log "[4/7] Cloning LAMMPS (Stable 2Aug2023_update3)..."
cd /tmp
if ! git clone --depth 1 --branch stable_2Aug2023_update3 https://github.com/lammps/lammps.git; then
    error "Git clone failed."
fi
cd lammps

#==============================================================================
# STEP 5: THE "ULTIMATE" COMPILATION (FIXED & ROBUST)
#==============================================================================
log "[5/7] Configuring CMake..."

# 1. Clean the battlefield
rm -rf build
mkdir -p build
cd build

# 2. Interrogate the MPI Fortran wrapper for the CORRECT link flags
# This prevents "cannot find -lmpifort" or "undefined reference" errors.
log "Detecting MPI Fortran Linker Flags..."
if command -v mpifort &> /dev/null; then
    MPI_FORT_LDFLAGS=$(mpifort --showme:link 2>/dev/null || mpifort -link_info 2>/dev/null || echo "-lmpi_mpifh")
    log "Detected Fortran Flags: ${MPI_FORT_LDFLAGS}"
else
    error "mpifort not found! Did you install libopenmpi-dev?"
fi

# 3. Export compiler wrappers
export CC=mpicc
export CXX=mpicxx
export FC=mpifort

cmake ../cmake -D CMAKE_BUILD_TYPE=Release -D CMAKE_CXX_COMPILER=mpicxx \
    -D CMAKE_CXX_FLAGS="-O3 -march=native -DNDEBUG" -D BUILD_MPI=yes -D BUILD_OMP=yes \
    -D PKG_OPENMP=yes -D PKG_INTEL=yes -D PKG_OPT=yes \
    -D PKG_MEAM=yes -D PKG_MANYBODY=yes -D PKG_PTM=yes -D PKG_VORONOI=yes -D DOWNLOAD_VORONOI=yes \
    -D PKG_REPLICA=yes -D PKG_MC=yes -D PKG_PYTHON=yes \
    -D PKG_PLUMED=yes -D DOWNLOAD_PLUMED=yes \
    -D PKG_IO=yes -D PKG_MPIIO=yes -D PKG_COMPRESS=yes \
    -D PKG_EXTRA-COMPUTE=yes -D PKG_EXTRA-DUMP=yes -D PKG_EXTRA-FIX=yes \
    -D WITH_GZIP=yes -D WITH_JPEG=yes -D WITH_PNG=yes -D WITH_FFMPEG=yes

# Note: The below is an alternative CMake command with """"ALL""""" packages enabled.

# cmake ../cmake -D CMAKE_BUILD_TYPE=Release -D CMAKE_INSTALL_PREFIX=/usr/local -D CMAKE_EXE_LINKER_FLAGS="${MPI_FORT_LDFLAGS}" -D BUILD_MPI=yes -D BUILD_OMP=yes -D BUILD_SHARED_LIBS=yes -D FFT=FFTW3 -D PKG_OPENMP=yes -D PKG_INTEL=yes -D PKG_OPT=yes -D PKG_KOKKOS=yes -D Kokkos_ENABLE_OPENMP=yes -D Kokkos_ENABLE_SERIAL=yes \
#     -D PKG_KSPACE=yes -D PKG_MANYBODY=yes -D PKG_MEAM=yes -D PKG_MOLECULE=yes -D PKG_RIGID=yes -D PKG_GRANULAR=yes -D PKG_ASPHERE=yes -D PKG_COLLOID=yes -D PKG_SNAP=yes -D PKG_REAXFF=yes -D PKG_QEQ=yes -D PKG_PERI=yes -D PKG_SPIN=yes -D PKG_DIPOLE=yes -D PKG_CLASS2=yes -D PKG_SHOCK=yes -D PKG_DIELECTRIC=yes \
#     -D PKG_KIM=yes -D DOWNLOAD_KIM=yes -D PKG_ML-SNAP=yes -D PKG_ML-PACE=yes -D PKG_ML-RANN=yes -D PKG_ML-IAP=yes -D PKG_ML-POD=yes -D PKG_MC=yes -D PKG_REPLICA=yes -D PKG_PYTHON=yes -D PKG_VORONOI=yes -D DOWNLOAD_VORONOI=yes -D PKG_PLUMED=yes -D DOWNLOAD_PLUMED=yes -D PKG_SCAFACOS=yes -D DOWNLOAD_SCAFACOS=yes \
#     -D PKG_COLVARS=yes -D PKG_PHONON=yes -D PKG_DIFFRACTION=yes -D PKG_PTM=yes -D PKG_IO=yes -D PKG_NETCDF=yes -D PKG_H5MD=yes -D PKG_MPIIO=yes -D PKG_COMPRESS=yes -D PKG_MISC=yes -D PKG_EXTRA-COMPUTE=yes -D PKG_EXTRA-DUMP=yes -D PKG_EXTRA-FIX=yes -D PKG_EXTRA-PAIR=yes \
#     -D WITH_GZIP=yes -D WITH_PNG=yes -D WITH_JPEG=yes -D WITH_FFMPEG=yes

log "Compiling with ${NPROC} cores..."
make -j"${NPROC}"

log "Installing binaries..."
sudo make install

#==============================================================================
# STEP 6: PYTHON INTERFACE & FINAL SETUP
#==============================================================================
log "[6/7] Installing Python Interface..."

# Ensure ldconfig knows about the new libs immediately
echo "/usr/local/lib" | sudo tee /etc/ld.so.conf.d/lammps.conf > /dev/null
sudo ldconfig

# Set variables so Python build can find the shared lib we just installed
export LD_LIBRARY_PATH=/usr/local/lib:${LD_LIBRARY_PATH:-}
export LAMMPS_LIB_PATH=/usr/local/lib

cd ../python
if ! pip3 install --break-system-packages . ; then
    warn "Python install failed. Attempting to install via make..."
    cd ../build
    sudo make install-python || error "Python interface installation failed."
fi

#==============================================================================
# STEP 7: BASH CONFIGURATION
#==============================================================================
log "[7/7] Configuring environment..."

# Idempotent Bashrc setup
if ! grep -q "LAMMPS_ULTIMATE" ~/.bashrc; then
cat >> ~/.bashrc << 'EOF'

# --- LAMMPS_ULTIMATE SETUP ---
export LD_LIBRARY_PATH=/usr/local/lib:$LD_LIBRARY_PATH
# Ensure Python finds the LAMMPS module
export PYTHONPATH=$PYTHONPATH:/usr/local/lib/python3/dist-packages

# Aliases
alias lmp='lmp'
alias lmp_mpi='mpirun -np $(nproc) lmp'
# Optimized execution for Intel CPUs (AVX + Mixed Precision)
alias lmp_intel='mpirun -np $(nproc) lmp -sf intel -pk intel 0 mode mixed'
EOF
fi

# Cleanup
cd ~
sudo rm -rf /tmp/lammps

echo ""
echo "========================================================"
echo -e "${GREEN}   LAMMPS INSTALLATION COMPLETE${NC}"
echo "========================================================"
echo " Architecture: Intel i7 (Optimized for AVX)"
echo " Version:      Stable 2Aug2023_update3"
echo " Python:       Enabled (import lammps)"
echo "========================================================"
echo " USAGE:"
echo -e "   ${YELLOW}source ~/.bashrc${NC}"
echo -e "   ${YELLOW}lmp_intel -in input.script${NC}"
echo ""
echo " NOTE: If you use potentials like EAM, standard MPI is fine."
echo " If you use ReaxFF or costly pairs, lmp_intel will fly."
echo "========================================================"