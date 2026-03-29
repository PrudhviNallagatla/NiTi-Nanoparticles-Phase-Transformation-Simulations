#!/usr/bin/env bash
set -euo pipefail

# Batch runner for Taguchi L9 inputs.
# Edit values below, save, then run: ./run_taguchi_l9.sh

PYTHON_BIN="python3"
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
INPUT_CSV="${SCRIPT_DIR}/data/taguchi_l9_inputs.csv"
OUTPUT_CSV="${SCRIPT_DIR}/taguchi_l9_predictions.csv"
MODEL_JSON=""

# Observation settings used for percentile range output.
SEM_COUNT=10
AGING_HOURS=0.0
AGGLOMERATION_RATE=0.08
GEO_SIGMA=1.35

export INPUT_CSV OUTPUT_CSV MODEL_JSON SEM_COUNT AGING_HOURS AGGLOMERATION_RATE GEO_SIGMA

"${PYTHON_BIN}" - <<'PY'
import csv
from pathlib import Path

from predictor import NitiNanoPredictor, ProcessParams

script_dir = Path.cwd()

input_csv = script_dir / "data" / "taguchi_l9_inputs.csv"
output_csv = script_dir / "taguchi_l9_predictions.csv"
model_json = ""
sem_count = 10
aging_hours = 0.0
agglomeration_rate = 0.08
geo_sigma = 1.35

# Read values from shell environment overrides if present.
import os
input_csv = Path(os.environ.get("INPUT_CSV", str(input_csv)))
output_csv = Path(os.environ.get("OUTPUT_CSV", str(output_csv)))
model_json = os.environ.get("MODEL_JSON", model_json)
sem_count = int(os.environ.get("SEM_COUNT", sem_count))
aging_hours = float(os.environ.get("AGING_HOURS", aging_hours))
agglomeration_rate = float(os.environ.get("AGGLOMERATION_RATE", agglomeration_rate))
geo_sigma = float(os.environ.get("GEO_SIGMA", geo_sigma))

predictor = NitiNanoPredictor.load(model_json) if model_json else NitiNanoPredictor()

rows_out = []
with input_csv.open("r", newline="", encoding="utf-8") as fh:
    reader = csv.DictReader(fh)
    for row in reader:
        proc = ProcessParams.from_experimental(
            voltage_v=float(row["voltage_v"]),
            current_a=float(row["current_a"]),
            duty_percent=float(row["duty_percent"]),
            frequency_khz=float(row["frequency_khz"]),
        )
        pred = predictor.predict_with_observation(
            proc=proc,
            sem_particle_count=sem_count,
            aging_hours=aging_hours,
            agglomeration_rate_per_hour=agglomeration_rate,
            geometric_sigma=geo_sigma,
        )
        rows_out.append(
            {
                "run": row["run"],
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

output_csv.parent.mkdir(parents=True, exist_ok=True)
with output_csv.open("w", newline="", encoding="utf-8") as fh:
    writer = csv.DictWriter(
        fh,
        fieldnames=[
            "run",
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
    writer.writerows(rows_out)

print(f"Wrote predicted table: {output_csv}")
PY

echo

echo "Done. Predicted table: ${OUTPUT_CSV}"
