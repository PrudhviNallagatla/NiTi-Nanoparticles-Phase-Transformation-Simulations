import glob
import re
import os
import sys
import numpy as np
import matplotlib
matplotlib.use('Agg') 
import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator, AutoMinorLocator, MaxNLocator
from scipy.signal import savgol_filter
from ovito.io import import_file
from ovito.modifiers import (
    PolyhedralTemplateMatchingModifier,
    ClusterAnalysisModifier,
    AtomicStrainModifier,
    ExpressionSelectionModifier
)

INPUT_DIR = "/home/rimuru/workspace/data/raw_output/niti_10nm_Ni50_single_crystal/thermal_cycle"
OUTPUT_DIR = INPUT_DIR.replace("raw_output", "analysis_output")

if not os.path.exists(INPUT_DIR):
    print(f"ERROR: Input directory does not exist!\n{INPUT_DIR}")
    sys.exit(1)

if not os.path.exists(OUTPUT_DIR):
    os.makedirs(OUTPUT_DIR)
    print(f"I: Created output directory: {OUTPUT_DIR}")

T_LOW  = 200.0
T_HIGH = 500.0
RATE   = 0.35
DT     = 0.002
RMSD   = 0.15
CORE_FRAC = 0.8

        
DELTA_T    = T_HIGH - T_LOW
RUN_RAMP   = round(DELTA_T / RATE / DT)
RUN_HOLD   = 50000
END_STAGE0 = RUN_HOLD
END_STAGE1 = END_STAGE0 + RUN_RAMP
END_STAGE2 = END_STAGE1 + RUN_HOLD
END_STAGE3 = END_STAGE2 + RUN_RAMP

def get_temp(step):
    if step <= END_STAGE0: return T_HIGH
    elif step <= END_STAGE1: return T_HIGH - ((step-END_STAGE0)/RUN_RAMP)*(T_HIGH-T_LOW)
    elif step <= END_STAGE2: return T_LOW
    elif step <= END_STAGE3: return T_LOW + ((step-END_STAGE2)/RUN_RAMP)*(T_HIGH-T_LOW)
    else: return T_HIGH

print(f"I: Analyzing data in: {INPUT_DIR}")
files = sorted(glob.glob(os.path.join(INPUT_DIR, "dump.cycle.*.gz")), 
               key=lambda f: int(re.search(r'dump\.cycle\.(\d+)', f).group(1)))

if not files:
    print("ERROR: No dump files found.")
    sys.exit(1)

                 
print("I: Configuring OVITO Modifiers...")

ptm = PolyhedralTemplateMatchingModifier(rmsd_cutoff=RMSD, output_orientation=False)

                                                
ref_pipeline = import_file(files[0]) 
strain_mod = AtomicStrainModifier(cutoff=3.5)
strain_mod.reference = ref_pipeline.source 

                                     
                                             
select_mod = ExpressionSelectionModifier(expression = "StructureType != 3")

                                
cluster_mod = ClusterAnalysisModifier(cutoff=3.5, sort_by_size=True, only_selected=True)

                
pipeline = import_file(files[0])
pipeline.modifiers.append(ptm)
pipeline.modifiers.append(strain_mod)
pipeline.modifiers.append(select_mod)
pipeline.modifiers.append(cluster_mod)

results = [] 

print(f"I: Processing {len(files)} frames (OVITO Engine)...")

for i, fname in enumerate(files):
    step = int(re.search(r'dump\.cycle\.(\d+)', fname).group(1))
    temp = get_temp(step)
    
                
    pipeline.source.load(fname)
    data = pipeline.compute()
    
              
    pos = data.particles.position
    com = np.mean(pos, axis=0)
    r_sq = np.sum((pos - com)**2, axis=1)
    max_r2 = np.max(r_sq)
    radius = np.sqrt(max_r2)
    is_core = r_sq < (CORE_FRAC * radius)**2
    is_surf = ~is_core
    
            
    structs = data.particles['Structure Type']
    n_total = len(structs)
    f_aus = np.count_nonzero(structs == 3) / n_total
    f_aus_c = np.count_nonzero(structs[is_core] == 3) / np.sum(is_core) if np.sum(is_core) > 0 else 0
    f_aus_s = np.count_nonzero(structs[is_surf] == 3) / np.sum(is_surf) if np.sum(is_surf) > 0 else 0
    
            
    if 'c_pe_all' in data.particles:
        pe_arr = data.particles['c_pe_all']
        pe_core = np.mean(pe_arr[is_core]) if np.sum(is_core) > 0 else 0
    else: pe_core = 0.0
        
               
    shear = data.particles['Shear Strain']
    shear_avg = np.mean(shear)
    
    stress_gpa = 0.0
    if 'c_st' in data.particles:
        st_tensor = data.particles['c_st']
        trace_sum = np.sum((st_tensor[:,0] + st_tensor[:,1] + st_tensor[:,2])[is_core])
        vol_core = (4/3) * np.pi * (radius * CORE_FRAC)**3
        stress_gpa = (trace_sum / 3.0) / vol_core / 10000.0
        
                    
    if 'ClusterAnalysis.cluster_count' in data.attributes:
        n_clusters = data.attributes['ClusterAnalysis.cluster_count']
    else:
        n_clusters = 0
    
               
    rg = np.sqrt(np.mean(r_sq))
    
    results.append([temp, f_aus, f_aus_c, f_aus_s, pe_core, shear_avg, stress_gpa, n_clusters, rg, step])
    
    if i % 50 == 0: print(f"   Step {step}: T={temp:.0f}K | Aus={f_aus:.2f}")

                
arr = np.array(results)
arr = arr[arr[:,9].argsort()]

               
np.savetxt(f"{OUTPUT_DIR}/extracted_data.txt", arr, header="Temp f_Aus f_Aus_Core f_Aus_Surf PE_Core Shear Stress Clusters Rg Step")
