from __future__ import annotations

import argparse
import csv
import json
import math
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Iterable, Optional

K_B = 1.380649e-23
N_A = 6.02214076e23


@dataclass
class ProcessParams:
    """Micro-EDM process inputs used by the predictor."""

    voltage_v: float
    capacitance_f: float
    duty_cycle: float
    pulse_on_us: float
    pulse_frequency_hz: float = 50_000.0
    current_a: Optional[float] = None

    @classmethod
    def from_experimental(
        cls,
        voltage_v: float,
        current_a: float,
        duty_percent: float,
        frequency_khz: float,
        capacitance_f: Optional[float] = None,
    ) -> "ProcessParams":
        """Build process parameters directly from typical Micro-EDM experiment columns."""
        duty_cycle = duty_percent / 100.0
        pulse_frequency_hz = frequency_khz * 1_000.0
        if pulse_frequency_hz <= 0:
            raise ValueError("frequency_khz must be > 0")
        if not (0 < duty_cycle < 1):
            raise ValueError("duty_percent must be between 0 and 100")

        period_s = 1.0 / pulse_frequency_hz
        pulse_on_s = duty_cycle * period_s
        pulse_on_us = pulse_on_s * 1e6

        # If capacitance is missing, infer an equivalent capacitance from pulse energy.
        # E_pulse = V * I * t_on, and E_pulse = 0.5 * C_eq * V^2 -> C_eq = 2 I t_on / V
        if capacitance_f is None:
            if voltage_v <= 0:
                raise ValueError("voltage_v must be > 0")
            capacitance_f = (2.0 * current_a * pulse_on_s) / voltage_v

        return cls(
            voltage_v=voltage_v,
            capacitance_f=capacitance_f,
            duty_cycle=duty_cycle,
            pulse_on_us=pulse_on_us,
            pulse_frequency_hz=pulse_frequency_hz,
            current_a=current_a,
        )

    def validate(self) -> None:
        if self.voltage_v <= 0:
            raise ValueError("voltage_v must be > 0")
        if self.capacitance_f <= 0:
            raise ValueError("capacitance_f must be > 0")
        if not (0 < self.duty_cycle < 1):
            raise ValueError("duty_cycle must be between 0 and 1")
        if self.pulse_on_us <= 0:
            raise ValueError("pulse_on_us must be > 0")
        if self.pulse_frequency_hz <= 0:
            raise ValueError("pulse_frequency_hz must be > 0")
        if self.current_a is not None and self.current_a <= 0:
            raise ValueError("current_a must be > 0 when provided")


@dataclass
class ModelParams:
    """Thermophysical and model constants (editable when calibrating)."""

    eta_anode: float = 0.094
    eta_additional_losses: float = 0.65
    h_vap_effective_j_per_kg: float = 8.0e6
    niti_density_liquid_kg_m3: float = 6050.0
    niti_molar_mass_kg_per_mol: float = 0.05328
    monomer_diameter_nm: float = 0.35
    dielectric_viscosity_pa_s: float = 0.02
    plasma_temperature_k: float = 10_000.0
    dielectric_temperature_k: float = 300.0
    freeze_temperature_k: float = 1100.0
    quench_tau_s: float = 2.5e-6
    quench_duty_factor: float = 2.0
    plasma_resistance_ohm: float = 0.02
    plasma_radius_coeff_m: float = 1.2e-5
    fractal_dimension: float = 1.8
    coagulation_scale: float = 1.0
    minimum_size_nm: float = 1.0


@dataclass
class Calibration:
    """Simple affine correction learned from experiments."""

    gain: float = 1.0
    bias_nm: float = 0.0


@dataclass
class PredictionResult:
    size_nm_raw: float
    size_nm_calibrated: float
    vapor_mass_kg: float
    coagulation_time_s: float
    bubble_radius_m: float
    number_density_m3: float
    size_nm_p10: float
    size_nm_p50: float
    size_nm_p90: float
    metadata: dict


class NitiNanoPredictor:
    """Physics-first Micro-EDM NiTi nanoparticle size predictor."""

    def __init__(
        self,
        model_params: Optional[ModelParams] = None,
        calibration: Optional[Calibration] = None,
    ) -> None:
        self.model_params = model_params or ModelParams()
        self.calibration = calibration or Calibration()

    def predict(self, proc: ProcessParams) -> PredictionResult:
        proc.validate()
        mp = self.model_params

        pulse_on_s = proc.pulse_on_us * 1e-6
        period_s = 1.0 / proc.pulse_frequency_hz
        pulse_off_s = max(0.0, period_s - pulse_on_s)

        e_cap_j = 0.5 * proc.capacitance_f * proc.voltage_v**2
        if proc.current_a is not None:
            # Experimental mode: energy per pulse comes from measured V-I-ton.
            e_pulse_j = proc.voltage_v * proc.current_a * pulse_on_s
            e_to_anode_j = e_pulse_j * mp.eta_anode * mp.eta_additional_losses
            energy_source = "vi_ton"
        else:
            # Legacy mode: use equivalent capacitor energy scaled by duty.
            e_to_anode_j = (
                e_cap_j * mp.eta_anode * mp.eta_additional_losses * proc.duty_cycle
            )
            energy_source = "capacitor"
        vapor_mass_kg = max(e_to_anode_j / mp.h_vap_effective_j_per_kg, 1e-21)

        current_a = (
            proc.current_a
            if proc.current_a is not None
            else proc.voltage_v / mp.plasma_resistance_ohm
        )
        bubble_radius_m = (
            mp.plasma_radius_coeff_m * (current_a**0.43) * (proc.pulse_on_us**0.44)
        )
        bubble_radius_m = max(bubble_radius_m, 2e-7)
        bubble_volume_m3 = (4.0 / 3.0) * math.pi * bubble_radius_m**3

        monomer_mass_kg = mp.niti_molar_mass_kg_per_mol / N_A
        n_particles = vapor_mass_kg / monomer_mass_kg
        number_density_m3 = max(n_particles / bubble_volume_m3, 1e6)

        quench_tau = mp.quench_tau_s * (1.0 + mp.quench_duty_factor * proc.duty_cycle)
        coagulation_time_s = min(
            self._time_to_temperature(proc, mp, quench_tau),
            pulse_off_s if pulse_off_s > 0 else period_s,
        )
        coagulation_time_s = max(coagulation_time_s, 1e-9)

        t_mid = 0.5 * coagulation_time_s
        t_mid_k = self._temperature_at(t_mid, mp, quench_tau)
        beta = (
            mp.coagulation_scale
            * (8.0 * K_B * t_mid_k)
            / (3.0 * mp.dielectric_viscosity_pa_s)
        )

        n_t = number_density_m3 / (1.0 + beta * number_density_m3 * coagulation_time_s)
        n_t = max(n_t, 1.0)

        mass_ratio = number_density_m3 / n_t
        primary_diameter_nm = mp.monomer_diameter_nm * mass_ratio ** (1.0 / 3.0)

        d_eff = 3.0 if t_mid_k >= mp.freeze_temperature_k else mp.fractal_dimension
        hydro_diameter_nm = primary_diameter_nm * mass_ratio ** (
            1.0 / d_eff - 1.0 / 3.0
        )
        size_nm_raw = max(hydro_diameter_nm, mp.minimum_size_nm)

        size_nm_cal = self.calibration.gain * size_nm_raw + self.calibration.bias_nm
        size_nm_cal = max(size_nm_cal, mp.minimum_size_nm)

        # Distribution-aware post-processing for SEM-sized datasets.
        # For small particle counts (for example, 5-15), sampling error can dominate.
        obs_stats = self._observation_adjusted_range(
            size_nm_base=size_nm_cal,
            sem_particle_count=10,
            aging_hours=0.0,
            agglomeration_rate_per_hour=0.08,
            geometric_sigma=1.35,
            proc=proc,
        )

        return PredictionResult(
            size_nm_raw=size_nm_raw,
            size_nm_calibrated=size_nm_cal,
            vapor_mass_kg=vapor_mass_kg,
            coagulation_time_s=coagulation_time_s,
            bubble_radius_m=bubble_radius_m,
            number_density_m3=number_density_m3,
            size_nm_p10=obs_stats["p10"],
            size_nm_p50=obs_stats["p50"],
            size_nm_p90=obs_stats["p90"],
            metadata={
                "energy_capacitor_j": e_cap_j,
                "energy_to_anode_j": e_to_anode_j,
                "energy_source": energy_source,
                "mid_temperature_k": t_mid_k,
                "mass_ratio": mass_ratio,
                "effective_fractal_dimension": d_eff,
                "aging_multiplier": obs_stats["aging_multiplier"],
                "sigma_ln_total": obs_stats["sigma_ln_total"],
                "sem_particle_count": obs_stats["sem_particle_count"],
            },
        )

    def predict_with_observation(
        self,
        proc: ProcessParams,
        sem_particle_count: int,
        aging_hours: float,
        agglomeration_rate_per_hour: float = 0.08,
        geometric_sigma: float = 1.35,
    ) -> PredictionResult:
        """Predict with explicit SEM sampling and sample-aging assumptions."""
        result = self.predict(proc)
        obs_stats = self._observation_adjusted_range(
            size_nm_base=result.size_nm_calibrated,
            sem_particle_count=sem_particle_count,
            aging_hours=aging_hours,
            agglomeration_rate_per_hour=agglomeration_rate_per_hour,
            geometric_sigma=geometric_sigma,
            proc=proc,
        )
        result.size_nm_p10 = obs_stats["p10"]
        result.size_nm_p50 = obs_stats["p50"]
        result.size_nm_p90 = obs_stats["p90"]
        result.metadata["aging_multiplier"] = obs_stats["aging_multiplier"]
        result.metadata["sigma_ln_total"] = obs_stats["sigma_ln_total"]
        result.metadata["sem_particle_count"] = obs_stats["sem_particle_count"]
        result.metadata["geometric_sigma"] = geometric_sigma
        result.metadata["agglomeration_rate_per_hour"] = agglomeration_rate_per_hour
        result.metadata["aging_hours"] = aging_hours
        return result

    @staticmethod
    def _observation_adjusted_range(
        size_nm_base: float,
        sem_particle_count: int,
        aging_hours: float,
        agglomeration_rate_per_hour: float,
        geometric_sigma: float,
        proc: ProcessParams,
    ) -> dict:
        sem_particle_count = max(int(sem_particle_count), 2)
        aging_hours = max(float(aging_hours), 0.0)
        geometric_sigma = max(float(geometric_sigma), 1.01)
        agglomeration_rate_per_hour = max(float(agglomeration_rate_per_hour), 0.0)

        duty_factor = 1.0 + 0.5 * proc.duty_cycle
        aging_multiplier = math.exp(
            agglomeration_rate_per_hour * duty_factor * aging_hours
        )
        p50 = max(size_nm_base * aging_multiplier, 1e-9)

        sigma_ln_intrinsic = math.log(geometric_sigma)
        sigma_ln_sampling = 0.6 / math.sqrt(sem_particle_count)
        sigma_ln_total = math.sqrt(sigma_ln_intrinsic**2 + sigma_ln_sampling**2)

        z10 = -1.2815515655446004
        z90 = 1.2815515655446004
        p10 = p50 * math.exp(z10 * sigma_ln_total)
        p90 = p50 * math.exp(z90 * sigma_ln_total)

        return {
            "p10": max(p10, 1e-9),
            "p50": p50,
            "p90": max(p90, 1e-9),
            "aging_multiplier": aging_multiplier,
            "sigma_ln_total": sigma_ln_total,
            "sem_particle_count": sem_particle_count,
        }

    def fit_calibration(self, records: Iterable[tuple[float, float]]) -> Calibration:
        """Fit affine correction using (predicted_raw_nm, measured_nm)."""
        x_vals = []
        y_vals = []
        for x_raw, y_meas in records:
            if x_raw > 0 and y_meas > 0:
                x_vals.append(float(x_raw))
                y_vals.append(float(y_meas))

        if len(x_vals) < 2:
            raise ValueError("Need at least two valid records for calibration")

        n = len(x_vals)
        x_mean = sum(x_vals) / n
        y_mean = sum(y_vals) / n
        s_xx = sum((x - x_mean) ** 2 for x in x_vals)
        if s_xx == 0:
            raise ValueError("Calibration failed: zero variance in predictor values")
        s_xy = sum((x - x_mean) * (y - y_mean) for x, y in zip(x_vals, y_vals))

        gain = s_xy / s_xx
        bias = y_mean - gain * x_mean
        self.calibration = Calibration(gain=gain, bias_nm=bias)
        return self.calibration

    def save(self, path: str | Path) -> None:
        payload = {
            "model_params": asdict(self.model_params),
            "calibration": asdict(self.calibration),
        }
        Path(path).write_text(json.dumps(payload, indent=2), encoding="utf-8")

    @classmethod
    def load(cls, path: str | Path) -> "NitiNanoPredictor":
        payload = json.loads(Path(path).read_text(encoding="utf-8"))
        return cls(
            model_params=ModelParams(**payload.get("model_params", {})),
            calibration=Calibration(**payload.get("calibration", {})),
        )

    @staticmethod
    def _temperature_at(t_s: float, mp: ModelParams, tau_s: float) -> float:
        return mp.dielectric_temperature_k + (
            mp.plasma_temperature_k - mp.dielectric_temperature_k
        ) * math.exp(-t_s / tau_s)

    def _time_to_temperature(
        self, proc: ProcessParams, mp: ModelParams, tau_s: float
    ) -> float:
        if mp.freeze_temperature_k <= mp.dielectric_temperature_k:
            return 1e-9
        ratio = (mp.freeze_temperature_k - mp.dielectric_temperature_k) / (
            mp.plasma_temperature_k - mp.dielectric_temperature_k
        )
        ratio = min(max(ratio, 1e-9), 0.999999)
        t_freeze = -tau_s * math.log(ratio)

        # Scale with pulse energy density so stronger pulses preserve coagulation window longer.
        pulse_on_s = proc.pulse_on_us * 1e-6
        period_s = 1.0 / proc.pulse_frequency_hz
        e_cap_j = 0.5 * proc.capacitance_f * proc.voltage_v**2
        energy_scale = max(
            0.25, min(4.0, 1.0 + 0.4 * math.log10(max(e_cap_j, 1e-12) / 1e-7 + 1.0))
        )
        return min(t_freeze * energy_scale, period_s)


def _read_calibration_csv(csv_path: Path) -> list[tuple[float, float]]:
    rows: list[tuple[float, float]] = []
    with csv_path.open("r", newline="", encoding="utf-8") as fh:
        reader = csv.DictReader(fh)
        field_set = set(reader.fieldnames or [])
        legacy_required = {
            "voltage_v",
            "capacitance_f",
            "duty_cycle",
            "pulse_on_us",
            "measured_size_nm",
        }
        exp_required = {
            "voltage_v",
            "current_a",
            "duty_percent",
            "frequency_khz",
            "measured_size_nm",
        }
        mode = ""
        if legacy_required.issubset(field_set):
            mode = "legacy"
        elif exp_required.issubset(field_set):
            mode = "experimental"
        else:
            raise ValueError(
                "Calibration CSV must contain either legacy columns "
                "(voltage_v, capacitance_f, duty_cycle, pulse_on_us, measured_size_nm) "
                "or experimental columns "
                "(voltage_v, current_a, duty_percent, frequency_khz, measured_size_nm)."
            )

        predictor = NitiNanoPredictor()
        for row in reader:
            if mode == "legacy":
                proc = ProcessParams(
                    voltage_v=float(row["voltage_v"]),
                    capacitance_f=float(row["capacitance_f"]),
                    duty_cycle=float(row["duty_cycle"]),
                    pulse_on_us=float(row["pulse_on_us"]),
                    pulse_frequency_hz=float(row.get("pulse_frequency_hz", 50_000.0)),
                )
            else:
                cap_raw = row.get("capacitance_f", "")
                cap_val = float(cap_raw) if cap_raw else None
                proc = ProcessParams.from_experimental(
                    voltage_v=float(row["voltage_v"]),
                    current_a=float(row["current_a"]),
                    duty_percent=float(row["duty_percent"]),
                    frequency_khz=float(row["frequency_khz"]),
                    capacitance_f=cap_val,
                )
            pred = predictor.predict(proc).size_nm_raw
            rows.append((pred, float(row["measured_size_nm"])))
    return rows


def _build_cli() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        description="Physics-first NiTi nanoparticle size predictor for Micro-EDM"
    )
    p.add_argument("--voltage", type=float, required=True, help="Open-gap voltage in V")
    p.add_argument(
        "--capacitance",
        type=float,
        default=0.0,
        help="Capacitance in F (legacy mode; optional in experimental mode)",
    )
    p.add_argument(
        "--duty", type=float, default=0.0, help="Duty cycle in (0, 1) for legacy mode"
    )
    p.add_argument(
        "--pulse-on-us",
        type=float,
        default=0.0,
        help="Pulse-on time in microseconds for legacy mode",
    )
    p.add_argument("--freq", type=float, default=50_000.0, help="Pulse frequency in Hz")
    p.add_argument(
        "--current",
        type=float,
        default=0.0,
        help="Measured current in A (experimental mode)",
    )
    p.add_argument(
        "--duty-percent",
        type=float,
        default=0.0,
        help="Duty cycle in percent (experimental mode)",
    )
    p.add_argument(
        "--freq-khz",
        type=float,
        default=0.0,
        help="Pulse frequency in kHz (experimental mode)",
    )
    p.add_argument("--model", type=str, default="", help="Load model JSON")
    p.add_argument(
        "--calibrate-csv", type=str, default="", help="CSV for affine calibration"
    )
    p.add_argument(
        "--sem-count",
        type=int,
        default=10,
        help="Number of nanoparticles identified in SEM for this run",
    )
    p.add_argument(
        "--aging-hours",
        type=float,
        default=0.0,
        help="Time delay between synthesis and SEM imaging in hours",
    )
    p.add_argument(
        "--agglomeration-rate",
        type=float,
        default=0.08,
        help="Effective agglomeration growth rate per hour",
    )
    p.add_argument(
        "--geo-sigma",
        type=float,
        default=1.35,
        help="Geometric sigma of underlying size distribution",
    )
    p.add_argument("--save-model", type=str, default="", help="Save model JSON")
    return p


def main() -> None:
    args = _build_cli().parse_args()
    predictor = (
        NitiNanoPredictor.load(args.model) if args.model else NitiNanoPredictor()
    )

    if args.calibrate_csv:
        cal_rows = _read_calibration_csv(Path(args.calibrate_csv))
        predictor.fit_calibration(cal_rows)

    experimental_mode = args.current > 0 and args.duty_percent > 0 and args.freq_khz > 0
    if experimental_mode:
        cap_val = args.capacitance if args.capacitance > 0 else None
        proc = ProcessParams.from_experimental(
            voltage_v=args.voltage,
            current_a=args.current,
            duty_percent=args.duty_percent,
            frequency_khz=args.freq_khz,
            capacitance_f=cap_val,
        )
    else:
        proc = ProcessParams(
            voltage_v=args.voltage,
            capacitance_f=args.capacitance,
            duty_cycle=args.duty,
            pulse_on_us=args.pulse_on_us,
            pulse_frequency_hz=args.freq,
        )
    result = predictor.predict_with_observation(
        proc=proc,
        sem_particle_count=args.sem_count,
        aging_hours=args.aging_hours,
        agglomeration_rate_per_hour=args.agglomeration_rate,
        geometric_sigma=args.geo_sigma,
    )

    payload = {
        "input_mode": "experimental" if experimental_mode else "legacy",
        "input": asdict(proc),
        "size_nm_raw": result.size_nm_raw,
        "size_nm_calibrated": result.size_nm_calibrated,
        "size_nm_p10": result.size_nm_p10,
        "size_nm_p50": result.size_nm_p50,
        "size_nm_p90": result.size_nm_p90,
        "vapor_mass_kg": result.vapor_mass_kg,
        "coagulation_time_s": result.coagulation_time_s,
        "bubble_radius_m": result.bubble_radius_m,
        "number_density_m3": result.number_density_m3,
        "model_params": asdict(predictor.model_params),
        "calibration": asdict(predictor.calibration),
        "metadata": result.metadata,
    }
    print(json.dumps(payload, indent=2))

    if args.save_model:
        predictor.save(args.save_model)


if __name__ == "__main__":
    main()
