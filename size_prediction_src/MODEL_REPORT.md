# NiTi Nanoparticle Size Prediction Model Report

## 1. Objective
This model predicts NiTi nanoparticle size from micro-EDM process parameters in a physics-informed way, then maps the prediction to an observation range suitable for SEM datasets with small sample counts and possible post-synthesis agglomeration.

Primary process inputs:
- Gap voltage $V$ (V)
- Current $I$ (A)
- Duty cycle $D$ (%)
- Frequency $f$ (kHz)

Primary outputs:
- Raw theoretical size (`size_nm_raw`)
- Calibrated size (`size_nm_calibrated`)
- Percentile prediction range (`size_nm_p10`, `size_nm_p50`, `size_nm_p90`)

## 2. Model Equations

### 2.1 Pulse timing
Convert frequency and duty cycle:

$$
 f_{Hz} = 1000 f_{kHz}, \quad T = \frac{1}{f_{Hz}}, \quad t_{on}=D_{frac}T, \quad D_{frac}=\frac{D}{100}
$$

### 2.2 Energy source term
The model supports two modes:

1. Experimental mode (preferred for your data):

$$
E_{pulse}=VIt_{on}
$$

2. Equivalent capacitor mode (legacy compatibility):

$$
E_{cap}=\frac{1}{2}CV^2
$$

Energy reaching the anode (effective deposition):

$$
E_{anode}=\eta_{anode}\eta_{loss}E_{source}
$$

### 2.3 Vaporized mass estimate

$$
m_{vap}=\frac{E_{anode}}{H_{vap,eff}}
$$

where $H_{vap,eff}$ is an effective vaporization enthalpy term that lumps latent + process losses.

### 2.4 Bubble and number density estimate
Bubble radius surrogate:

$$
R_b = k_r I^{0.43} t_{on,\mu s}^{0.44}
$$

Bubble volume:

$$
V_b = \frac{4}{3}\pi R_b^3
$$

Monomer mass:

$$
m_{mono}=\frac{M_{NiTi}}{N_A}
$$

Initial particle number density:

$$
N_0=\frac{m_{vap}/m_{mono}}{V_b}
$$

### 2.5 Coagulation window and cooling
Cooling time constant:

$$
\tau_q = \tau_0(1+k_D D_{frac})
$$

Temperature decay model:

$$
T(t)=T_{amb} + (T_{plasma}-T_{amb})e^{-t/\tau_q}
$$

Characteristic coagulation duration is capped by pulse timing and freeze threshold behavior in code.

### 2.6 Brownian coagulation closure
Collision kernel surrogate:

$$
\beta = s_c\frac{8k_BT}{3\mu}
$$

Population reduction:

$$
N_t = \frac{N_0}{1+\beta N_0 t_c}
$$

Mass growth ratio:

$$
\Gamma = \frac{N_0}{N_t}
$$

Primary diameter:

$$
d_p = d_{mono}\Gamma^{1/3}
$$

Hydrodynamic growth with effective fractal dimension $d_f$:

$$
d_h = d_p\Gamma^{(1/d_f - 1/3)}
$$

This gives the raw model size `size_nm_raw`.

### 2.7 Affine calibration
The calibrated size is:

$$
y_{cal}=g\,y_{raw}+b
$$

where $g$ and $b$ are fitted from data if calibration is used.

### 2.8 Observation model (SEM uncertainty + aging/agglomeration)
Aging multiplier:

$$
A_{age}=\exp(r_{agg}(1+0.5D_{frac})t_{age,h})
$$

Median prediction:

$$
P50=A_{age}\,y_{cal}
$$

Log-normal spread:

$$
\sigma_{ln,intr} = \ln(\sigma_g), \quad
\sigma_{ln,samp} = \frac{0.6}{\sqrt{n_{SEM}}}, \quad
\sigma_{ln} = \sqrt{\sigma_{ln,intr}^2+\sigma_{ln,samp}^2}
$$

Percentiles:

$$
P10=P50\exp(-1.28155\sigma_{ln}), \quad
P90=P50\exp(+1.28155\sigma_{ln})
$$

## 3. Files and Folder Layout
- `predictor.py`: core mathematical model and CLI
- `run_predictor.sh`: one-click runner with editable parameters
- `data/taguchi_l9_inputs.csv`: process design matrix
- `data/experimental_responses.csv`: SEM + zeta response table
- `data/optimized_validation.csv`: optimized run validation record
- `refernce.md`: extracted source literature notes

## 4. Step-by-Step Run Guide

### Step 1: Open the shell script
Edit parameter block in `run_predictor.sh`:
- `VOLTAGE_V`
- `CURRENT_A`
- `DUTY_PERCENT`
- `FREQ_KHZ`
- and optional SEM/agglomeration parameters.

### Step 2: Make executable (first time only)
```bash
chmod +x run_predictor.sh
```

### Step 3: Run prediction
```bash
./run_predictor.sh
```

### Step 4: Read output
The script writes a JSON output to:
- `prediction_output.json`

This contains central and percentile predictions ready for reporting.

## 5. Publication-Style Workflow Recommendation
1. Use `data/taguchi_l9_inputs.csv` as benchmark process conditions.
2. Generate predictions for each run with fixed model constants.
3. Report model assumptions explicitly (Section 2 equations).
4. Compare predicted central tendency to SEM average and D50.
5. Compare predicted P10/P90 to experimental D10/D90.
6. Report MAE, RMSE, and MAPE in a results table.
7. Use ML only as a secondary enhancement after this mechanistic baseline.

## 6. Scope and Limits
- This is a reduced-order physics-informed model, not a full CFD/plasma solver.
- Parameters like effective vaporization enthalpy and kernel scale are lumped terms.
- Small SEM counts and delayed imaging can dominate uncertainty; therefore percentile outputs are preferred over single-number claims.

## 7. Next Technical Upgrade Path
1. Fit calibration and agglomeration parameters against your own measurements.
2. Add uncertainty propagation (Monte Carlo over uncertain constants).
3. Add residual ML correction model with cross-validation.
4. Validate on out-of-sample experiments beyond the L9 grid.
