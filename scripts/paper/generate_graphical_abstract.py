#!/usr/bin/env python3
"""
Generate a graphical abstract using actual OVITO renderings plus flowchart overlays.

Primary output is SVG (Inkscape editable). A PNG preview is also generated.
"""

import argparse
import os

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
from matplotlib.patches import Circle, FancyArrowPatch, FancyBboxPatch

REPO_ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))
OUT_DIR = os.path.join(REPO_ROOT, "paper", "figures")

ASSETS = {
    "5nm_start": os.path.join(OUT_DIR, "5nm_se_start.png"),
    "5nm_mid": os.path.join(OUT_DIR, "5nm_se_middle.png"),
    "5nm_end": os.path.join(OUT_DIR, "5nm_se_finish.png"),
    "10nm_start": os.path.join(OUT_DIR, "10nm_se_start.png"),
    "10nm_mid": os.path.join(OUT_DIR, "10nm_se_middle.png"),
    "10nm_end": os.path.join(OUT_DIR, "10nm_se_finish.png"),
}

PALETTE = {
    "bg": "#f4f7fb",
    "panel": "#ffffff",
    "edge": "#8fa2b8",
    "text": "#16202a",
    "muted": "#44546a",
    "arrow": "#1663a6",
    "accent": "#0f3b66",
    "ok": "#2f855a",
    "warn": "#c05621",
}


def ensure_assets():
    missing = [p for p in ASSETS.values() if not os.path.exists(p)]
    if missing:
        msg = "\n".join(missing)
        raise FileNotFoundError(f"Missing OVITO assets:\n{msg}")


def style_image_axis(ax, title=None):
    ax.set_xticks([])
    ax.set_yticks([])
    for spine in ax.spines.values():
        spine.set_visible(True)
        spine.set_linewidth(1.2)
        spine.set_color(PALETTE["edge"])
    if title:
        ax.set_title(
            title, fontsize=10.5, color=PALETTE["accent"], pad=6, fontweight="bold"
        )


def add_row_flow_arrow(canvas_ax, x0, x1, y):
    canvas_ax.add_patch(
        FancyArrowPatch(
            (x0, y),
            (x1, y),
            arrowstyle="-|>",
            mutation_scale=13,
            lw=2.0,
            color=PALETTE["arrow"],
            zorder=6,
        )
    )


def draw_flowchart_boxes(canvas_ax):
    box_thermal = FancyBboxPatch(
        (0.05, 0.835),
        0.24,
        0.09,
        boxstyle="round,pad=0.008,rounding_size=0.02",
        linewidth=1.1,
        edgecolor=PALETTE["edge"],
        facecolor="#edf5ff",
    )
    box_mech = FancyBboxPatch(
        (0.33, 0.835),
        0.30,
        0.09,
        boxstyle="round,pad=0.008,rounding_size=0.02",
        linewidth=1.1,
        edgecolor=PALETTE["edge"],
        facecolor="#edf5ff",
    )
    box_size = FancyBboxPatch(
        (0.68, 0.835),
        0.27,
        0.09,
        boxstyle="round,pad=0.008,rounding_size=0.02",
        linewidth=1.1,
        edgecolor=PALETTE["edge"],
        facecolor="#edf5ff",
    )
    canvas_ax.add_patch(box_thermal)
    canvas_ax.add_patch(box_mech)
    canvas_ax.add_patch(box_size)

    canvas_ax.text(
        0.17,
        0.88,
        "Thermal cycling\n200-500 K",
        ha="center",
        va="center",
        fontsize=10.2,
        color=PALETTE["accent"],
        fontweight="bold",
    )
    canvas_ax.text(
        0.48,
        0.88,
        "Mechanical loading at 500 K\n(SE path from OVITO)",
        ha="center",
        va="center",
        fontsize=10.2,
        color=PALETTE["accent"],
        fontweight="bold",
    )
    canvas_ax.text(
        0.815,
        0.88,
        "Size effect in spatial response\ncore-surface decoupling",
        ha="center",
        va="center",
        fontsize=10.2,
        color=PALETTE["accent"],
        fontweight="bold",
    )

    canvas_ax.add_patch(
        FancyArrowPatch(
            (0.29, 0.88),
            (0.33, 0.88),
            arrowstyle="-|>",
            mutation_scale=14,
            lw=2.0,
            color=PALETTE["arrow"],
        )
    )
    canvas_ax.add_patch(
        FancyArrowPatch(
            (0.63, 0.88),
            (0.68, 0.88),
            arrowstyle="-|>",
            mutation_scale=14,
            lw=2.0,
            color=PALETTE["arrow"],
        )
    )


def add_takeaway_and_callout(canvas_ax):
    callout = FancyBboxPatch(
        (0.80, 0.29),
        0.18,
        0.23,
        boxstyle="round,pad=0.01,rounding_size=0.02",
        linewidth=1.1,
        edgecolor=PALETTE["edge"],
        facecolor="#fffaf0",
    )
    canvas_ax.add_patch(callout)
    canvas_ax.text(
        0.89,
        0.47,
        "10 nm: transformed core",
        ha="center",
        va="center",
        fontsize=9.5,
        color=PALETTE["warn"],
        fontweight="bold",
    )
    canvas_ax.text(
        0.89,
        0.42,
        "+ retained-B2 surface rim",
        ha="center",
        va="center",
        fontsize=9.5,
        color=PALETTE["ok"],
        fontweight="bold",
    )
    canvas_ax.text(
        0.89,
        0.37,
        "surface exclusion zone",
        ha="center",
        va="center",
        fontsize=8.6,
        color=PALETTE["muted"],
    )
    canvas_ax.text(
        0.89,
        0.33,
        "(~1-2 nm)",
        ha="center",
        va="center",
        fontsize=8.6,
        color=PALETTE["muted"],
    )

    canvas_ax.text(
        0.5,
        0.08,
        "Thermal pathway is suppressed at nanoscale, but stress-driven reversible transformation remains active, with stronger core-surface decoupling at 10 nm.",
        ha="center",
        va="center",
        fontsize=11,
        color=PALETTE["text"],
        fontweight="bold",
    )


def build_figure(width=16, height=9):
    ensure_assets()

    fig = plt.figure(figsize=(width, height), facecolor=PALETTE["bg"])
    gs = GridSpec(
        nrows=8,
        ncols=6,
        figure=fig,
        left=0.04,
        right=0.78,
        top=0.80,
        bottom=0.16,
        wspace=0.05,
        hspace=0.14,
    )

    canvas_ax = fig.add_axes([0, 0, 1, 1])
    canvas_ax.set_xlim(0, 1)
    canvas_ax.set_ylim(0, 1)
    canvas_ax.axis("off")

    canvas_ax.text(
        0.5,
        0.965,
        "Graphical Abstract: Nanoscale NiTi shows thermal suppression but stress-driven reversibility",
        ha="center",
        va="top",
        fontsize=16.2,
        color=PALETTE["text"],
        fontweight="bold",
    )

    draw_flowchart_boxes(canvas_ax)

    # 5 nm row
    ax_5_start = fig.add_subplot(gs[0:4, 0:2])
    ax_5_mid = fig.add_subplot(gs[0:4, 2:4])
    ax_5_end = fig.add_subplot(gs[0:4, 4:6])

    img_5_start = plt.imread(ASSETS["5nm_start"])
    img_5_mid = plt.imread(ASSETS["5nm_mid"])
    img_5_end = plt.imread(ASSETS["5nm_end"])

    ax_5_start.imshow(img_5_start)
    ax_5_mid.imshow(img_5_mid)
    ax_5_end.imshow(img_5_end)

    style_image_axis(ax_5_start, "5 nm: start")
    style_image_axis(ax_5_mid, "5 nm: peak strain")
    style_image_axis(ax_5_end, "5 nm: recovered")

    # 10 nm row
    ax_10_start = fig.add_subplot(gs[4:8, 0:2])
    ax_10_mid = fig.add_subplot(gs[4:8, 2:4])
    ax_10_end = fig.add_subplot(gs[4:8, 4:6])

    img_10_start = plt.imread(ASSETS["10nm_start"])
    img_10_mid = plt.imread(ASSETS["10nm_mid"])
    img_10_end = plt.imread(ASSETS["10nm_end"])

    ax_10_start.imshow(img_10_start)
    ax_10_mid.imshow(img_10_mid)
    ax_10_end.imshow(img_10_end)

    style_image_axis(ax_10_start, "10 nm: start")
    style_image_axis(ax_10_mid, "10 nm: peak strain")
    style_image_axis(ax_10_end, "10 nm: recovered")

    # Row labels
    canvas_ax.text(
        0.022,
        0.58,
        "5 nm",
        rotation=90,
        ha="center",
        va="center",
        fontsize=12,
        color=PALETTE["accent"],
        fontweight="bold",
    )
    canvas_ax.text(
        0.022,
        0.315,
        "10 nm",
        rotation=90,
        ha="center",
        va="center",
        fontsize=12,
        color=PALETTE["accent"],
        fontweight="bold",
    )

    # Flow arrows along each row
    add_row_flow_arrow(canvas_ax, 0.245, 0.31, 0.57)
    add_row_flow_arrow(canvas_ax, 0.485, 0.55, 0.57)
    add_row_flow_arrow(canvas_ax, 0.245, 0.31, 0.305)
    add_row_flow_arrow(canvas_ax, 0.485, 0.55, 0.305)

    # Highlight zone in 10nm recovered state
    ring = Circle(
        (0.67, 0.24), 0.035, fill=False, lw=2.0, ec=PALETTE["warn"], linestyle="--"
    )
    canvas_ax.add_patch(ring)
    canvas_ax.add_patch(
        FancyArrowPatch(
            (0.705, 0.24),
            (0.81, 0.36),
            arrowstyle="-|>",
            mutation_scale=13,
            lw=1.9,
            color=PALETTE["warn"],
        )
    )

    add_takeaway_and_callout(canvas_ax)
    return fig


def parse_args():
    parser = argparse.ArgumentParser(
        description="Generate OVITO-based graphical abstract as SVG and PNG"
    )
    parser.add_argument(
        "--svg",
        default=os.path.join(OUT_DIR, "graphical_abstract.svg"),
        help="Path for SVG output",
    )
    parser.add_argument(
        "--png",
        default=os.path.join(OUT_DIR, "graphical_abstract.png"),
        help="Path for PNG preview output",
    )
    parser.add_argument(
        "--dpi",
        type=int,
        default=300,
        help="DPI for PNG output",
    )
    return parser.parse_args()


def main():
    args = parse_args()
    os.makedirs(os.path.dirname(args.svg), exist_ok=True)
    os.makedirs(os.path.dirname(args.png), exist_ok=True)

    fig = build_figure()
    fig.savefig(args.svg, format="svg", bbox_inches="tight", pad_inches=0.03)
    fig.savefig(
        args.png, format="png", dpi=args.dpi, bbox_inches="tight", pad_inches=0.03
    )
    plt.close(fig)

    print(f"[OK] Wrote SVG: {args.svg}")
    print(f"[OK] Wrote PNG: {args.png}")


if __name__ == "__main__":
    main()
