import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from scipy import integrate, stats
from scipy.ndimage import gaussian_filter1d
from scipy.signal import find_peaks, savgol_filter
import os

# Set Times New Roman as default font for all matplotlib plots
plt.rcParams["font.family"] = "Times New Roman"

# --- PARAMETERS ---
# Windows file path - use raw string (r"") or forward slashes
FILE_PATH = (
    r"D:\MD Sims\SME_SE\data\analysis_output\uncategorized\5nm_1grain_hysteresis.txt"
)

# Alternative ways to specify the path (choose one):
# FILE_PATH = "d:/MD Sims/SME_SE/data/raw_output/niti_2nm_Ni50_amorphous/5nm_1grain_hysteresis.txt"
# FILE_PATH = os.path.join("d:", "MD Sims", "SME_SE", "data", "raw_output", "niti_2nm_Ni50_amorphous", "5nm_1grain_hysteresis.txt")

SMOOTHING_WINDOW = 50
GAUSSIAN_SIGMA = 5
DERIVATIVE_THRESHOLD = 0.01  # For detecting phase transformation onset

# Verify file exists
if not os.path.exists(FILE_PATH):
    print(f"❌ ERROR: File not found at: {FILE_PATH}")
    print(f"Current working directory: {os.getcwd()}")
    print("\nPlease check:")
    print("1. File path is correct")
    print("2. File name spelling")
    print("3. Directory exists")
    exit()
else:
    print(f"✓ File found: {FILE_PATH}")

# 1. Load the Data (.txt file with space-separated values)
try:
    # Read the .txt file - it has space-separated values
    df = pd.read_csv(
        FILE_PATH,
        sep=r"\s+",  # Matches any whitespace (spaces/tabs)
        comment="#",  # Skip lines starting with #
        header=None,  # No header row in the file
        skipinitialspace=True,  # Skip spaces after delimiter
    )

    # Assign column names based on your file structure
    df.columns = ["step", "T_global", "T_core", "PE", "Vol", "Austenite_Frac"]
    print(f"✓ Loaded {len(df)} data points")
    print(
        f"✓ Temperature range: {df['T_global'].min():.1f} - {df['T_global'].max():.1f} K"
    )
    print(
        f"✓ Austenite fraction range: {df['Austenite_Frac'].min():.4f} - {df['Austenite_Frac'].max():.4f}"
    )

except FileNotFoundError:
    print(f"❌ ERROR: Could not find file: {FILE_PATH}")
    exit()
except Exception as e:
    print(f"❌ ERROR reading file: {e}")
    print("\nFirst few lines of file (for debugging):")
    try:
        with open(FILE_PATH, "r") as f:
            for i, line in enumerate(f):
                if i < 5:
                    print(line.strip())
                else:
                    break
    except:
        pass
    exit()

# 2. Smooth the Data
df["T_smooth"] = df["T_global"].rolling(window=SMOOTHING_WINDOW, center=True).mean()
df["PE_smooth"] = df["PE"].rolling(window=SMOOTHING_WINDOW, center=True).mean()
df["Vol_smooth"] = df["Vol"].rolling(window=SMOOTHING_WINDOW, center=True).mean()
df["Aust_smooth"] = (
    df["Austenite_Frac"].rolling(window=SMOOTHING_WINDOW, center=True).mean()
)
df = df.dropna()

print(f"✓ After smoothing: {len(df)} valid points")

# 3. Separate Cooling and Heating Curves
min_temp_idx = df["T_smooth"].idxmin()
cooling_df = df.loc[:min_temp_idx].copy().sort_values("T_smooth", ascending=False)
heating_df = df.loc[min_temp_idx:].copy().sort_values("T_smooth", ascending=True)

print(f"✓ Cooling phase: {len(cooling_df)} points")
print(f"✓ Heating phase: {len(heating_df)} points")

# ========== TRANSFORMATION TEMPERATURE EXTRACTION ==========


def find_transformation_temps(temp, austenite_frac, phase_name="Martensite"):
    """
    Extract Ms, Mf (cooling) or As, Af (heating) using multiple methods

    Methods:
    1. Tangent Method (ASTM Standard)
    2. Maximum Derivative Method
    3. 10%-90% Fraction Method
    4. Inflection Point Detection
    """

    # Ensure monotonic temperature for interpolation
    if len(temp) < 10:
        return None, None, None

    # Convert to numpy arrays for easier handling
    temp_arr = temp.values if hasattr(temp, "values") else np.array(temp)
    frac_arr = (
        austenite_frac.values
        if hasattr(austenite_frac, "values")
        else np.array(austenite_frac)
    )

    # Calculate derivative of austenite fraction w.r.t temperature
    if len(frac_arr) >= 11:
        dAust_dT = np.gradient(
            savgol_filter(frac_arr, window_length=11, polyorder=3), temp_arr
        )
    else:
        dAust_dT = np.gradient(frac_arr, temp_arr)

    # METHOD 1: Maximum Derivative (Peak dA/dT)
    if phase_name == "Martensite":
        # For cooling: looking for maximum decrease (minimum dA/dT)
        idx_max_deriv = np.argmin(dAust_dT)
    else:
        # For heating: looking for maximum increase (maximum dA/dT)
        idx_max_deriv = np.argmax(dAust_dT)

    T_peak = temp_arr[idx_max_deriv]

    # METHOD 2: 10%-90% Fraction Method (ASTM E2004)
    T_start = None
    T_finish = None

    try:
        if phase_name == "Martensite":
            # Cooling: Ms at 10% martensite (90% austenite), Mf at 90% martensite (10% austenite)
            # Need to check if values reach these thresholds
            if frac_arr.min() <= 0.10 and frac_arr.max() >= 0.90:
                T_start = np.interp(0.90, frac_arr[::-1], temp_arr[::-1])
                T_finish = np.interp(0.10, frac_arr[::-1], temp_arr[::-1])
        else:
            # Heating: As at 10% austenite, Af at 90% austenite
            if frac_arr.min() <= 0.10 and frac_arr.max() >= 0.90:
                T_start = np.interp(0.10, frac_arr, temp_arr)
                T_finish = np.interp(0.90, frac_arr, temp_arr)
    except Exception as e:
        print(f"⚠ Warning in {phase_name} interpolation: {e}")

    # METHOD 3: Tangent Method (requires clear inflection)
    T_tangent = None
    try:
        window = min(20, len(temp_arr) // 4)  # Adaptive window size
        max_slope = 0

        for i in range(window, len(temp_arr) - window):
            local_slope = abs(
                np.polyfit(
                    temp_arr[i - window : i + window],
                    frac_arr[i - window : i + window],
                    1,
                )[0]
            )
            if local_slope > max_slope:
                max_slope = local_slope
                T_tangent = temp_arr[i]
    except:
        pass

    return {
        "T_peak_derivative": T_peak,
        "T_start_10pct": T_start,
        "T_finish_90pct": T_finish,
        "T_tangent": T_tangent,
        "max_derivative_value": dAust_dT[idx_max_deriv],
    }


# Extract Transformation Temperatures
print("\n" + "=" * 70)
print("Analyzing transformation temperatures...")
print("=" * 70)

cooling_results = find_transformation_temps(
    cooling_df["T_smooth"], cooling_df["Aust_smooth"], "Martensite"
)

heating_results = find_transformation_temps(
    heating_df["T_smooth"], heating_df["Aust_smooth"], "Austenite"
)

# Assign standard nomenclature
Ms = cooling_results["T_start_10pct"]
Mf = cooling_results["T_finish_90pct"]
As = heating_results["T_start_10pct"]
Af = heating_results["T_finish_90pct"]

# Alternative estimates
Ms_alt = cooling_results["T_peak_derivative"]
Mf_alt = cooling_results["T_tangent"]
As_alt = heating_results["T_peak_derivative"]
Af_alt = heating_results["T_tangent"]

# ========== ADDITIONAL PROPERTIES ==========


# Heat Capacity Derivative (to detect latent heat peak)
def calculate_heat_capacity(temp, pe):
    pe_smooth = gaussian_filter1d(pe, sigma=GAUSSIAN_SIGMA)
    dPE_dT = np.gradient(pe_smooth, temp)
    return dPE_dT


cp_cooling = calculate_heat_capacity(
    cooling_df["T_smooth"].values, cooling_df["PE_smooth"].values
)
cp_heating = calculate_heat_capacity(
    heating_df["T_smooth"].values, heating_df["PE_smooth"].values
)

# Hysteresis Width
if Ms and As:
    hysteresis_width = abs(Ms - As)
else:
    hysteresis_width = None

# Transformation Range
if Ms and Mf:
    martensite_range = abs(Ms - Mf)
else:
    martensite_range = None

if As and Af:
    austenite_range = abs(Af - As)
else:
    austenite_range = None

# Energy Dissipation
T_common = np.linspace(df["T_smooth"].min(), df["T_smooth"].max(), 500)
PE_cool_interp = np.interp(T_common, cooling_df["T_smooth"], cooling_df["PE_smooth"])
PE_heat_interp = np.interp(
    T_common, heating_df["T_smooth"][::-1], heating_df["PE_smooth"][::-1]
)
hysteresis_area = integrate.trapz(PE_heat_interp - PE_cool_interp, T_common)

# ========== VISUALIZATION ==========

fig, axes = plt.subplots(2, 3, figsize=(20, 12))
fig.suptitle(
    "NiTi 5nm Transformation Temperature Analysis",
    fontsize=16,
    fontweight="bold",
    family="Times New Roman",
)

# Plot 1: Austenite Fraction vs Temperature (PRIMARY ANALYSIS)
axes[0, 0].plot(
    cooling_df["T_smooth"],
    cooling_df["Aust_smooth"],
    "b-",
    linewidth=2.5,
    label="Cooling",
)
axes[0, 0].plot(
    heating_df["T_smooth"],
    heating_df["Aust_smooth"],
    "r-",
    linewidth=2.5,
    label="Heating",
)

# Mark transformation temperatures
if Ms:
    axes[0, 0].axvline(
        Ms, color="blue", linestyle="--", linewidth=2, label=f"Ms = {Ms:.1f} K"
    )
    axes[0, 0].plot(Ms, 0.90, "bo", markersize=10)
if Mf:
    axes[0, 0].axvline(
        Mf, color="cyan", linestyle="--", linewidth=2, label=f"Mf = {Mf:.1f} K"
    )
    axes[0, 0].plot(Mf, 0.10, "co", markersize=10)
if As:
    axes[0, 0].axvline(
        As, color="red", linestyle="--", linewidth=2, label=f"As = {As:.1f} K"
    )
    axes[0, 0].plot(As, 0.10, "ro", markersize=10)
if Af:
    axes[0, 0].axvline(
        Af, color="orange", linestyle="--", linewidth=2, label=f"Af = {Af:.1f} K"
    )
    axes[0, 0].plot(Af, 0.90, color="orange", marker="o", markersize=10)

axes[0, 0].set_xlabel(
    "Temperature (K)", fontsize=13, fontweight="bold", family="Times New Roman"
)
axes[0, 0].set_ylabel(
    "Austenite Fraction", fontsize=13, fontweight="bold", family="Times New Roman"
)
axes[0, 0].set_title(
    "Phase Transformation Hysteresis",
    fontweight="bold",
    fontsize=14,
    family="Times New Roman",
)
axes[0, 0].legend(loc="best", fontsize=10)
axes[0, 0].grid(True, alpha=0.3)
axes[0, 0].set_ylim(0, 1)

# Plot 2: Derivative of Austenite Fraction (Peak Detection)
dAust_dT_cool = np.gradient(cooling_df["Aust_smooth"], cooling_df["T_smooth"])
dAust_dT_heat = np.gradient(heating_df["Aust_smooth"], heating_df["T_smooth"])

axes[0, 1].plot(
    cooling_df["T_smooth"], dAust_dT_cool, "b-", linewidth=2, label="Cooling (dA/dT)"
)
axes[0, 1].plot(
    heating_df["T_smooth"], dAust_dT_heat, "r-", linewidth=2, label="Heating (dA/dT)"
)
axes[0, 1].axhline(0, color="black", linestyle="-", linewidth=0.5)

if Ms_alt:
    axes[0, 1].axvline(
        Ms_alt, color="blue", linestyle=":", label=f"Ms (peak) = {Ms_alt:.1f} K"
    )
if As_alt:
    axes[0, 1].axvline(
        As_alt, color="red", linestyle=":", label=f"As (peak) = {As_alt:.1f} K"
    )

axes[0, 1].set_xlabel(
    "Temperature (K)", fontsize=13, fontweight="bold", family="Times New Roman"
)
axes[0, 1].set_ylabel(
    "dA/dT (K⁻¹)", fontsize=13, fontweight="bold", family="Times New Roman"
)
axes[0, 1].set_title(
    "Transformation Rate Analysis",
    fontweight="bold",
    fontsize=14,
    family="Times New Roman",
)
axes[0, 1].legend(fontsize=10)
axes[0, 1].grid(True, alpha=0.3)

# Plot 3: Heat Capacity (Latent Heat Detection)
axes[0, 2].plot(
    cooling_df["T_smooth"].values[:-1],
    cp_cooling[:-1],
    "b-",
    linewidth=2,
    label="Cooling",
)
axes[0, 2].plot(
    heating_df["T_smooth"].values[:-1],
    cp_heating[:-1],
    "r-",
    linewidth=2,
    label="Heating",
)

if Ms:
    axes[0, 2].axvline(Ms, color="blue", linestyle="--", alpha=0.5)
if Mf:
    axes[0, 2].axvline(Mf, color="cyan", linestyle="--", alpha=0.5)
if As:
    axes[0, 2].axvline(As, color="red", linestyle="--", alpha=0.5)
if Af:
    axes[0, 2].axvline(Af, color="orange", linestyle="--", alpha=0.5)

axes[0, 2].set_xlabel(
    "Temperature (K)", fontsize=13, fontweight="bold", family="Times New Roman"
)
axes[0, 2].set_ylabel(
    "dPE/dT (eV/K)", fontsize=13, fontweight="bold", family="Times New Roman"
)
axes[0, 2].set_title(
    "Heat Capacity (Latent Heat)",
    fontweight="bold",
    fontsize=14,
    family="Times New Roman",
)
axes[0, 2].legend()
axes[0, 2].grid(True, alpha=0.3)

# Plot 4: PE vs T (Energy Hysteresis)
axes[1, 0].plot(
    cooling_df["T_smooth"], cooling_df["PE_smooth"], "b-", linewidth=2, label="Cooling"
)
axes[1, 0].plot(
    heating_df["T_smooth"], heating_df["PE_smooth"], "r-", linewidth=2, label="Heating"
)
axes[1, 0].fill_between(
    T_common,
    PE_cool_interp,
    PE_heat_interp,
    alpha=0.3,
    color="gray",
    label=f"Hysteresis = {hysteresis_area:.2f} eV·K",
)
axes[1, 0].set_xlabel(
    "Temperature (K)", fontsize=13, fontweight="bold", family="Times New Roman"
)
axes[1, 0].set_ylabel(
    "Potential Energy (eV)", fontsize=13, fontweight="bold", family="Times New Roman"
)
axes[1, 0].set_title(
    "Energy Hysteresis Loop", fontweight="bold", fontsize=14, family="Times New Roman"
)
axes[1, 0].legend()
axes[1, 0].grid(True, alpha=0.3)

# Plot 5: Volume vs Temperature
axes[1, 1].plot(
    cooling_df["T_smooth"], cooling_df["Vol_smooth"], "b-", linewidth=2, label="Cooling"
)
axes[1, 1].plot(
    heating_df["T_smooth"], heating_df["Vol_smooth"], "r-", linewidth=2, label="Heating"
)
axes[1, 1].set_xlabel(
    "Temperature (K)", fontsize=13, fontweight="bold", family="Times New Roman"
)
axes[1, 1].set_ylabel(
    "Volume (Å³)", fontsize=13, fontweight="bold", family="Times New Roman"
)
axes[1, 1].set_title(
    "Thermal Expansion", fontweight="bold", fontsize=14, family="Times New Roman"
)
axes[1, 1].legend()
axes[1, 1].grid(True, alpha=0.3)

# Plot 6: Transformation Summary Diagram
axes[1, 2].axis("off")
summary_text = f"""
Transformation Temperatures (K)

Martensite Start (Ms):  {f'{Ms:.1f}' if Ms else 'N/A':>8}
Martensite Finish (Mf): {f'{Mf:.1f}' if Mf else 'N/A':>8}

Austenite Start (As):   {f'{As:.1f}' if As else 'N/A':>8}
Austenite Finish (Af):  {f'{Af:.1f}' if Af else 'N/A':>8}


Transformation Characteristics

Hysteresis Width:       {f'{hysteresis_width:.1f}' if hysteresis_width else 'N/A':>8} K
Martensite Range:       {f'{martensite_range:.1f}' if martensite_range else 'N/A':>8} K
Austenite Range:        {f'{austenite_range:.1f}' if austenite_range else 'N/A':>8} K

Energy Dissipation:     {hysteresis_area:>8.2f} eV·K


Phase Transformation Status
"""

if Ms and Mf and As and Af:
    summary_text += "✓ Full Transformation Detected\n"
elif Ms or As:
    summary_text += "⚠ Partial Transformation\n"
else:
    summary_text += "❌ No Transformation (Size-Effect)\n"

axes[1, 2].text(
    0.1,
    0.5,
    summary_text,
    fontsize=11,
    family="Times New Roman",
    verticalalignment="center",
    fontweight="normal",
)

plt.tight_layout()

# Save in the same directory as the input file
output_dir = os.path.dirname(FILE_PATH)
output_png = os.path.join(output_dir, "transformation_temperatures.png")
plt.savefig(output_png, dpi=300, bbox_inches="tight")
print(f"\n✓ Generated: {output_png}")

# ========== DETAILED REPORT ==========

print("\n" + "=" * 70)
print("           Transformation Temperature Analysis Report           ")
print("=" * 70)

print("\n" + "=" * 70)
print("              Transformation Temperatures (ASTM Method)           ")
print("=" * 70)
print(
    f"  Ms (Martensite Start - 10% method):    {f'{Ms:.1f}' if Ms else 'Not Detected':>18} K"
)
print(
    f"  Mf (Martensite Finish - 90% method):   {f'{Mf:.1f}' if Mf else 'Not Detected':>18} K"
)
print(
    f"  As (Austenite Start - 10% method):     {f'{As:.1f}' if As else 'Not Detected':>18} K"
)
print(
    f"  Af (Austenite Finish - 90% method):    {f'{Af:.1f}' if Af else 'Not Detected':>18} K"
)
print("=" * 70)

print("\n" + "=" * 70)
print("            Alternative Detection Methods (Verification)          ")
print("=" * 70)
print(
    f"  Ms (Peak Derivative):                  {f'{Ms_alt:.1f}' if Ms_alt else 'N/A':>18} K"
)
print(
    f"  As (Peak Derivative):                  {f'{As_alt:.1f}' if As_alt else 'N/A':>18} K"
)
print(
    f"  Ms (Tangent Method):                   {f'{Mf_alt:.1f}' if Mf_alt else 'N/A':>18} K"
)
print(
    f"  As (Tangent Method):                   {f'{Af_alt:.1f}' if Af_alt else 'N/A':>18} K"
)
print("=" * 70)

print("\n" + "=" * 70)
print("                  Transformation Characteristics                  ")
print("=" * 70)
if hysteresis_width:
    print(f"  Thermal Hysteresis (Ms - As):          {hysteresis_width:>18.2f} K")
else:
    print(f"  Thermal Hysteresis:                    {'Not Calculable':>18}")

if martensite_range:
    print(f"  Martensite Formation Range (Ms-Mf):    {martensite_range:>18.2f} K")
else:
    print(f"  Martensite Range:                      {'Not Calculable':>18}")

if austenite_range:
    print(f"  Austenite Formation Range (Af-As):     {austenite_range:>18.2f} K")
else:
    print(f"  Austenite Range:                       {'Not Calculable':>18}")

print(f"  Energy Dissipation (Hysteresis Area):  {hysteresis_area:>18.2f} eV·K")
print("=" * 70)

print("\n" + "=" * 70)
print("                         Interpretation                            ")
print("=" * 70)

if Ms and Mf and As and Af:
    print("  ✓ Complete phase transformation detected")
    print("  ✓ System exhibits reversible martensitic transformation")
    print(
        f"  ✓ Transformation completeness: ~{(df['Austenite_Frac'].max() - df['Austenite_Frac'].min())*100:.1f}%"
    )
elif Ms or As:
    print("  ⚠ Partial phase transformation detected")
    print("  ⚠ Possible size-induced suppression of full transformation")
    print("  → Consider: Nucleation barriers in nanoparticles")
else:
    print("  ❌ No Transformation Detected")
    print("  → Size-Effect: Particle too small for phase transformation")
    print("  → Critical size for NiTi transformation: ~10-20 nm")
    print("  → Recommend: Test larger particles or different stoichiometry")

print("=" * 70)

# Save results to CSV
output_csv1 = os.path.join(output_dir, "transformation_temperatures.csv")
output_csv2 = os.path.join(output_dir, "transformation_characteristics.csv")

results_df = pd.DataFrame(
    {
        "Temperature": ["Ms", "Mf", "As", "Af", "Ms_alt", "As_alt"],
        "Value_K": [
            Ms if Ms else np.nan,
            Mf if Mf else np.nan,
            As if As else np.nan,
            Af if Af else np.nan,
            Ms_alt if Ms_alt else np.nan,
            As_alt if As_alt else np.nan,
        ],
        "Method": [
            "10% Austenite",
            "90% Martensite",
            "10% Austenite",
            "90% Austenite",
            "Peak Derivative",
            "Peak Derivative",
        ],
    }
)

characteristics_df = pd.DataFrame(
    {
        "Property": [
            "Hysteresis_Width",
            "Martensite_Range",
            "Austenite_Range",
            "Energy_Dissipation",
        ],
        "Value": [
            hysteresis_width if hysteresis_width else np.nan,
            martensite_range if martensite_range else np.nan,
            austenite_range if austenite_range else np.nan,
            hysteresis_area,
        ],
        "Unit": ["K", "K", "K", "eV·K"],
    }
)

results_df.to_csv(output_csv1, index=False)
characteristics_df.to_csv(output_csv2, index=False)

print(f"\n✓ Saved: {output_csv1}")
print(f"✓ Saved: {output_csv2}")
print("\n" + "=" * 70)
print("Analysis complete! Check the output directory for results.")
print("=" * 70)
