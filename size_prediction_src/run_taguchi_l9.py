#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
import json
import re
from dataclasses import asdict
from pathlib import Path
from typing import Any

from predictor import NitiNanoPredictor, ProcessParams


def _sanitize_label(value: str) -> str:
    label = re.sub(r"[^A-Za-z0-9_.-]+", "_", value.strip())
    label = label.strip("._-")
    return label or "run"


def _build_cli() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        description="Batch predictor for Taguchi L9 inputs with per-run JSON outputs"
    )
    script_dir = Path(__file__).resolve().parent
    p.add_argument(
        "--input-csv",
        type=Path,
        default=script_dir / "data" / "taguchi_l9_inputs.csv",
        help="Input CSV with columns: run, voltage_v, current_a, duty_percent, frequency_khz",
    )
    p.add_argument(
        "--output-dir",
        type=Path,
        default=script_dir / "output" / "taguchi_l9",
        help="Directory where per-run JSON files are written",
    )
    p.add_argument(
        "--output-csv",
        type=Path,
        default=script_dir / "output" / "taguchi_l9" / "taguchi_l9_predictions.csv",
        help="Summary CSV file path",
    )
    p.add_argument("--model", type=Path, default=None, help="Optional model JSON path")
    p.add_argument("--sem-count", type=int, default=10)
    p.add_argument("--aging-hours", type=float, default=0.0)
    p.add_argument("--agglomeration-rate", type=float, default=0.08)
    p.add_argument("--geo-sigma", type=float, default=1.35)
    return p


def _read_rows(input_csv: Path) -> list[dict[str, str]]:
    with input_csv.open("r", newline="", encoding="utf-8") as fh:
        reader = csv.DictReader(fh)
        required = {"run", "voltage_v", "current_a", "duty_percent", "frequency_khz"}
        headers = set(reader.fieldnames or [])
        missing = sorted(required - headers)
        if missing:
            raise ValueError(
                f"Input CSV missing required columns: {', '.join(missing)}"
            )
        return list(reader)


def _build_run_payload(
    row: dict[str, str],
    proc: ProcessParams,
    pred: Any,
    sem_count: int,
    aging_hours: float,
    agglomeration_rate: float,
    geo_sigma: float,
) -> dict[str, Any]:
    return {
        "run": row["run"],
        "input_mode": "experimental",
        "input": {
            "voltage_v": float(row["voltage_v"]),
            "current_a": float(row["current_a"]),
            "duty_percent": float(row["duty_percent"]),
            "frequency_khz": float(row["frequency_khz"]),
            "resolved_process_params": asdict(proc),
        },
        "prediction": {
            "size_nm_raw": pred.size_nm_raw,
            "size_nm_calibrated": pred.size_nm_calibrated,
            "size_nm_p10": pred.size_nm_p10,
            "size_nm_p50": pred.size_nm_p50,
            "size_nm_p90": pred.size_nm_p90,
            "vapor_mass_kg": pred.vapor_mass_kg,
            "coagulation_time_s": pred.coagulation_time_s,
            "bubble_radius_m": pred.bubble_radius_m,
            "number_density_m3": pred.number_density_m3,
        },
        "observation_settings": {
            "sem_count": sem_count,
            "aging_hours": aging_hours,
            "agglomeration_rate_per_hour": agglomeration_rate,
            "geo_sigma": geo_sigma,
        },
        "metadata": pred.metadata,
    }


def main() -> None:
    args = _build_cli().parse_args()

    input_csv = args.input_csv.resolve()
    output_dir = args.output_dir.resolve()
    output_csv = args.output_csv.resolve()

    predictor = (
        NitiNanoPredictor.load(args.model) if args.model else NitiNanoPredictor()
    )

    rows = _read_rows(input_csv)
    output_dir.mkdir(parents=True, exist_ok=True)
    output_csv.parent.mkdir(parents=True, exist_ok=True)

    summary_rows: list[dict[str, str]] = []
    for row in rows:
        run_label = _sanitize_label(row["run"])
        proc = ProcessParams.from_experimental(
            voltage_v=float(row["voltage_v"]),
            current_a=float(row["current_a"]),
            duty_percent=float(row["duty_percent"]),
            frequency_khz=float(row["frequency_khz"]),
        )
        pred = predictor.predict_with_observation(
            proc=proc,
            sem_particle_count=args.sem_count,
            aging_hours=args.aging_hours,
            agglomeration_rate_per_hour=args.agglomeration_rate,
            geometric_sigma=args.geo_sigma,
        )

        run_payload = _build_run_payload(
            row=row,
            proc=proc,
            pred=pred,
            sem_count=args.sem_count,
            aging_hours=args.aging_hours,
            agglomeration_rate=args.agglomeration_rate,
            geo_sigma=args.geo_sigma,
        )
        run_json_path = output_dir / f"{run_label}.json"
        run_json_path.write_text(json.dumps(run_payload, indent=2), encoding="utf-8")

        summary_rows.append(
            {
                "run": row["run"],
                "json_file": run_json_path.name,
                "voltage_v": row["voltage_v"],
                "current_a": row["current_a"],
                "duty_percent": row["duty_percent"],
                "frequency_khz": row["frequency_khz"],
                "pred_size_raw_nm": f"{pred.size_nm_raw:.6f}",
                "pred_size_calibrated_nm": f"{pred.size_nm_calibrated:.6f}",
                "pred_p10_nm": f"{pred.size_nm_p10:.6f}",
                "pred_p50_nm": f"{pred.size_nm_p50:.6f}",
                "pred_p90_nm": f"{pred.size_nm_p90:.6f}",
            }
        )

    with output_csv.open("w", newline="", encoding="utf-8") as fh:
        writer = csv.DictWriter(
            fh,
            fieldnames=[
                "run",
                "json_file",
                "voltage_v",
                "current_a",
                "duty_percent",
                "frequency_khz",
                "pred_size_raw_nm",
                "pred_size_calibrated_nm",
                "pred_p10_nm",
                "pred_p50_nm",
                "pred_p90_nm",
            ],
        )
        writer.writeheader()
        writer.writerows(summary_rows)

    print(f"Wrote summary table: {output_csv}")
    print(f"Wrote {len(summary_rows)} per-run JSON files in: {output_dir}")


if __name__ == "__main__":
    main()
