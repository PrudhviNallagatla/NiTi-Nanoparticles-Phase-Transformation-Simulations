#!/bin/bash
# ==========================================================================
# NiTi NANOPARTICLE STRUCTURE GENERATOR
# Purpose: Create physically realistic NiTi nanoparticles with defects
# Methods: ASE (single crystal) + Atomsk (polycrystalline grain boundaries)
# Output: LAMMPS data file + XYZ for visualization
# ==========================================================================

set -euo pipefail  # Exit on error, undefined variables, or pipe failures

# ============================================
# SCIENTIFIC PARAMETERS
# Purpose: Control particle geometry, composition, and microstructure
# ============================================

PROJECT_ROOT="/home/rimuru/workspace"
SEED=42  # Random seed for reproducible defect placement

# --- GEOMETRY ---
SIZE=10              # Particle diameter in nanometers
SHAPE="sphere"       # Options: sphere, ellipsoid, faceted, rough, blob

# --- COMPOSITION ---
NI_PERCENT=50.0      # Ni atomic percentage (rest is Ti)
                     # Note: Equiatomic (Ni50Ti50) is most studied composition

# --- COMPUTATIONAL RESOURCES ---
CORES=-1             # CPU cores for parallel processing (-1 = use all available)

# ============================================
# MICROSTRUCTURE CONTROL (CRITICAL)
# ============================================

# Grain Structure:
# - GRAIN_SIZE = 0    → Single crystal (no grain boundaries)
# - GRAIN_SIZE > 0    → Polycrystalline with specified grain size
# 
# Scientific Context:
#   • Grain boundaries affect transformation temperatures
#   • Smaller grains = more interfaces = altered properties
#   • Typical experimental range: 5-50 nm grains
GRAIN_SIZE=0.0       # Target grain diameter in nm (0 = single crystal)

# --- DEFECT CONCENTRATIONS ---
# Purpose: Introduce realistic imperfections found in real materials
VACANCIES=0.01           # Bulk vacancy concentration (1% = high defect density)
SURFACE_VACANCIES=0.03   # Surface vacancy concentration (3% = rough surface)
ANTISITES=0.005          # Antisite defects: Ni on Ti site or vice versa (0.5%)

# Scientific Notes:
#   • Vacancies enable atomic mobility during phase transformation
#   • Surface defects are nucleation sites for transformation
#   • Antisites create local compositional disorder

# ============================================
# TOOL VALIDATION (MANDATORY PRE-FLIGHT CHECK)
# Purpose: Ensure all required scientific software is installed
# ============================================

echo "========================================"
echo "VALIDATING SCIENTIFIC SOFTWARE"
echo "========================================"

# ----------------
# Check 1: Atomsk
# ----------------
# Purpose: Industry-standard tool for creating polycrystalline structures
# Citation: Hirel (2015) Comput. Phys. Commun. 197:212-219
if ! command -v atomsk &> /dev/null; then
    echo "❌ CRITICAL ERROR: Atomsk not found!"
    echo ""
    echo "Atomsk is REQUIRED for polycrystalline structures."
    echo "It is the peer-reviewed standard for grain boundary generation."
    echo ""
    echo "INSTALLATION (Linux):"
    echo "  1. Download: https://atomsk.univ-lille.fr/dl.php"
    echo "  2. Install:"
    echo "     wget https://atomsk.univ-lille.fr/code/atomsk_b0.13.1_Linux-x86-64.tar.gz"
    echo "     tar -xzf atomsk_b0.13.1_Linux-x86-64.tar.gz"
    echo "     sudo mv atomsk /usr/local/bin/"
    echo "     sudo chmod +x /usr/local/bin/atomsk"
    echo ""
    echo "CITATION (include in publications):"
    echo "  Hirel, P. (2015). Atomsk: A tool for manipulating and"
    echo "  converting atomic data files. Computer Physics Communications,"
    echo "  197, 212-219. DOI: 10.1016/j.cpc.2015.07.012"
    echo "========================================"
    exit 1
fi

# Display Atomsk version for logging
ATOMSK_VERSION=$(atomsk --version 2>&1 | head -n 1 || echo "unknown")
echo "✅ Atomsk found: ${ATOMSK_VERSION}"

# ----------------
# Check 2: Python Dependencies
# ----------------
# Required packages:
#   • ase:      Atomic Simulation Environment (structure manipulation)
#   • numpy:    Numerical arrays
#   • scipy:    Scientific computing
#   • psutil:   System resource monitoring
if ! python3 -c "import ase, numpy, scipy, psutil" 2>/dev/null; then
    echo "❌ ERROR: Missing Python dependencies"
    echo ""
    echo "Install with:"
    echo "  pip install ase numpy scipy tqdm threadpoolctl psutil"
    echo ""
    echo "CITATIONS:"
    echo "  • ASE: Larsen et al. (2017) J. Phys.: Condens. Matter 29:273002"
    echo "========================================"
    exit 1
fi

echo "✅ Python dependencies validated (ase, numpy, scipy, psutil)"
echo "========================================"

# ============================================
# STRUCTURE GENERATION SETUP
# ============================================

# Determine structure type based on grain size
if (( $(python3 -c "print(int(${GRAIN_SIZE} > 0))") )); then
    TYPE="nanocrystalline_${GRAIN_SIZE}nm_grains_atomsk"
    STRUCTURE_TYPE="polycrystalline"
    echo "Structure type: POLYCRYSTALLINE"
    echo "  Method: Atomsk Voronoi tessellation"
    echo "  Grain size: ${GRAIN_SIZE} nm"
else
    TYPE="single_crystal"
    STRUCTURE_TYPE="single_crystal"
    echo "Structure type: SINGLE CRYSTAL"
    echo "  Method: ASE direct lattice construction"
fi

# Generate output filenames
NI_INT=$(python3 -c "print(int(${NI_PERCENT}))")  # Integer composition for filename
TI_INT=$((100 - NI_INT))

OUTPUT="${PROJECT_ROOT}/data/structures/niti_${SIZE}nm_Ni${NI_INT}_${TYPE}.data"
OUTPUT_XYZ="${OUTPUT%.data}.xyz"

mkdir -p "${PROJECT_ROOT}/data/structures"

# ============================================
# DISPLAY GENERATION PARAMETERS
# ============================================
echo ""
echo "========================================"
echo "GENERATING NiTi NANOPARTICLE"
echo "========================================"
echo "Geometry:"
echo "  ├─ Diameter:   ${SIZE} nm"
echo "  └─ Shape:      ${SHAPE}"
echo ""
echo "Composition:"
echo "  └─ Formula:    Ni${NI_INT}Ti${TI_INT} (${NI_PERCENT}% Ni)"
echo ""
echo "Microstructure:"
echo "  ├─ Type:       ${STRUCTURE_TYPE}"

if [[ "${STRUCTURE_TYPE}" == "polycrystalline" ]]; then
    echo "  ├─ Grain size: ${GRAIN_SIZE} nm"
    echo "  ├─ Method:     Atomsk Voronoi tessellation"
    echo "  └─ Citation:   Hirel (2015) DOI: 10.1016/j.cpc.2015.07.012"
else
    echo "  └─ Method:     ASE B2 lattice (a=3.015 Å)"
fi

echo ""
echo "Defects:"
echo "  ├─ Bulk vacancies:     ${VACANCIES} ($(python3 -c "print(${VACANCIES}*100)")%)"
echo "  ├─ Surface vacancies:  ${SURFACE_VACANCIES} ($(python3 -c "print(${SURFACE_VACANCIES}*100)")%)"
echo "  └─ Antisite defects:   ${ANTISITES} ($(python3 -c "print(${ANTISITES}*100)")%)"
echo ""
echo "Random seed: ${SEED} (for reproducibility)"
echo "========================================"

# ============================================
# BUILD COMMAND LINE
# Purpose: Construct Python script call with all parameters
# ============================================

CMD="python3 \"${PROJECT_ROOT}/src/ase_np_gen.py\" \
    --diameter ${SIZE} \
    --shape ${SHAPE} \
    --ni-percent ${NI_PERCENT} \
    --output \"${OUTPUT}\" \
    --seed ${SEED} \
    --cores ${CORES} \
    --xyz"

# Add grain structure parameters ONLY for polycrystalline
if [[ "${STRUCTURE_TYPE}" == "polycrystalline" ]]; then
    CMD="${CMD} --grain-size ${GRAIN_SIZE}"
    CMD="${CMD} --use-atomsk"  # MANDATORY flag for grain boundary creation
fi

# Add defect parameters
CMD="${CMD} --vacancies ${VACANCIES}"
CMD="${CMD} --surface-vacancies ${SURFACE_VACANCIES}"
CMD="${CMD} --antisites ${ANTISITES}"

echo ""
echo "Executing Python generator..."
echo "Command: ${CMD}"
echo "========================================"
echo ""

# ============================================
# EXECUTE STRUCTURE GENERATION
# ============================================
eval ${CMD}

# ============================================
# POST-GENERATION VALIDATION
# Purpose: Verify output quality and correctness
# ============================================

# Check if Atomsk was actually invoked (for polycrystalline structures)
if [[ "${STRUCTURE_TYPE}" == "polycrystalline" ]]; then
    if [[ -f "${OUTPUT}" ]]; then
        # Look for Atomsk signature in output file comments
        if ! grep -q "Atomsk" "${OUTPUT}" 2>/dev/null; then
            echo ""
            echo "⚠️  WARNING: Atomsk signature not found in output!"
            echo "   This may indicate Atomsk failed silently."
            echo "   Review the log output above for error messages."
            echo ""
        else
            echo ""
            echo "✅ ATOMSK VALIDATION PASSED"
            echo "   Grain boundaries successfully generated"
        fi
    fi
fi

# ============================================
# SUCCESS SUMMARY
# ============================================
if [[ -f "${OUTPUT}" ]]; then
    echo ""
    echo "========================================"
    echo "✅ GENERATION SUCCESSFUL"
    echo "========================================"
    echo "Output files:"
    echo "  ├─ LAMMPS data: ${OUTPUT}"
    echo "  └─ XYZ format:  ${OUTPUT_XYZ}"
    echo ""
    
    # Display file sizes
    LAMMPS_SIZE=$(stat -c%s "${OUTPUT}" 2>/dev/null || stat -f%z "${OUTPUT}" 2>/dev/null || echo "0")
    XYZ_SIZE=$(stat -c%s "${OUTPUT_XYZ}" 2>/dev/null || stat -f%z "${OUTPUT_XYZ}" 2>/dev/null || echo "0")
    
    echo "File sizes:"
    echo "  ├─ LAMMPS: $(python3 -c "print(f'{${LAMMPS_SIZE}/1024:.1f} KB')")"
    echo "  └─ XYZ:    $(python3 -c "print(f'{${XYZ_SIZE}/1024:.1f} KB')")"
    echo ""
    
    # Display required citations
    echo "REQUIRED CITATIONS (include in publications):"
    if [[ "${STRUCTURE_TYPE}" == "polycrystalline" ]]; then
        echo "  ✅ Atomsk: Hirel (2015) Computer Physics Communications"
        echo "     197:212-219. DOI: 10.1016/j.cpc.2015.07.012"
    fi
    echo "  ✅ ASE: Larsen et al. (2017) Journal of Physics:"
    echo "     Condensed Matter 29:273002. DOI: 10.1088/1361-648X/aa680e"
    echo "========================================"
    echo ""
    
    # ============================================
    # NEXT STEPS & VALIDATION GUIDE
    # ============================================
    echo "NEXT STEPS:"
    echo ""
    echo "1. Visual inspection:"
    echo "   ovito ${OUTPUT_XYZ}"
    echo "   → Check particle shape and size"
    echo "   → Look for obvious artifacts"
    echo ""
    
    if [[ "${STRUCTURE_TYPE}" == "polycrystalline" ]]; then
        # Expected grain count (rough estimate: volume ratio)
        EXPECTED_GRAINS=$(python3 -c "print(max(1, int((${SIZE}/${GRAIN_SIZE})**3)))")
        
        echo "2. Validate grain structure (in OVITO):"
        echo "   a) Add modifier: Analysis → Grain Segmentation"
        echo "   b) Expected grain count: ~${EXPECTED_GRAINS} grains"
        echo "   c) Check grain boundary thickness: 2-4 Å is realistic"
        echo "   d) Verify grain size distribution is reasonable"
        echo ""
        echo "   If grain structure looks wrong:"
        echo "   → Check Atomsk warnings in log above"
        echo "   → Try different GRAIN_SIZE value"
        echo "   → Ensure GRAIN_SIZE << SIZE (needs at least 2-3 grains)"
        echo ""
    fi
    
    echo "3. Energy minimization (MANDATORY before use):"
    echo "   bash scripts/minimize_structure.sh"
    echo "   → Relaxes high-energy interfaces"
    echo "   → Required for stable MD simulations"
    echo ""
    
    echo "4. Structural validation:"
    echo "   → Count atoms: should be ~thousands for 10 nm particle"
    echo "   → Check composition: Ni/Ti ratio should match ${NI_INT}:${TI_INT}"
    echo "   → Verify defect density: use OVITO coordination analysis"
    echo ""
    echo "========================================"
    
else
    # Generation failed
    echo ""
    echo "❌ ERROR: Structure generation failed!"
    echo "========================================"
    echo "Troubleshooting:"
    echo "  1. Check Python script errors above"
    echo "  2. Verify ASE installation: python3 -c 'import ase'"
    if [[ "${STRUCTURE_TYPE}" == "polycrystalline" ]]; then
        echo "  3. Verify Atomsk is working: atomsk --version"
        echo "  4. Try with GRAIN_SIZE=0 (single crystal) to test ASE"
    fi
    echo "  5. Check write permissions in ${PROJECT_ROOT}/data/structures/"
    echo "========================================"
    exit 1
fi