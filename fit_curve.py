#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Curve parameter estimator for:
x = X + t*cos(theta) - exp(M*|t|)*sin(0.3*t)*sin(theta)
y = 42 + t*sin(theta) + exp(M*|t|)*sin(0.3*t)*cos(theta)

We only use numpy and the standard library.
Strategy:
1) Load CSV of (x,y) points.
2) Define a loss: rotate points by -theta after shifting by X
   => first coord is t_est, second coord is s_est
   => s_est should equal exp(M*|t|)*sin(0.3*t)
   Minimize median L1 residual + small penalty if t outside (6,60).
3) Do a coarse grid search (fast).
4) Do coordinate refinement around the coarse best (no SciPy).
5) Print final parameters and write a Desmos string to desmos_string.txt
"""

import sys
import math
import csv
from statistics import median
from typing import List, Tuple
import numpy as np
import matplotlib.pyplot as plt

# ----------------------------
# Helper: read xy csv (x,y per row; headers optional)
# ----------------------------
def load_xy(csv_path: str) -> np.ndarray:
    xs, ys = [], []
    with open(csv_path, "r", newline="") as f:
        reader = csv.reader(f)
        for row in reader:
            if not row:
                continue
            # Be tolerant of headers
            try:
                x = float(row[0])
                y = float(row[1])
            except ValueError:
                # skip header or bad line
                continue
            xs.append(x)
            ys.append(y)
    if not xs:
        raise RuntimeError("No numeric rows found in CSV.")
    return np.vstack([np.array(xs, dtype=float), np.array(ys, dtype=float)]).T  # shape (N,2)

# ----------------------------
# Core math for a batch of points
# ----------------------------
def rotate_and_split(points: np.ndarray, theta: float, X: float) -> Tuple[np.ndarray, np.ndarray]:
    """
    Shift by (X,42), then rotate by -theta.
    Returns (t_est, s_est) arrays of shape (N,).
    """
    c = math.cos(theta)
    s = math.sin(theta)
    u = points[:, 0] - X      # x - X
    v = points[:, 1] - 42.0   # y - 42

    # R(-theta) * [u; v]  => [ t_est ; s_est ]
    t_est =  c * u + s * v
    s_est = -s * u + c * v
    return t_est, s_est

def s_model(t: np.ndarray, M: float) -> np.ndarray:
    return np.exp(M * np.abs(t)) * np.sin(0.3 * t)

def loss(points: np.ndarray, theta: float, M: float, X: float, lam: float = 0.1) -> float:
    """
    L = median |s_est - s_model(t_est)|  + lam * boundary_penalty(t)
    """
    t_est, s_est = rotate_and_split(points, theta, X)
    s_hat = s_model(t_est, M)
    resid = np.abs(s_est - s_hat)
    # median L1 is robust to outliers and doesn't push to overfit
    base = float(median(resid.tolist()))

    # penalty if t is outside (6, 60)
    low = np.maximum(0.0, 6.0 - t_est)
    high = np.maximum(0.0, t_est - 60.0)
    penalty = float(np.mean(low + high))

    return base + lam * penalty

# ----------------------------
# Coarse search (very fast)
# ----------------------------
def coarse_search(points: np.ndarray) -> Tuple[float, float, float, float]:
    best = (1e9, None, None, None)  # (loss, theta, M, X)
    # theta in radians: 0..50 deg
    for theta_deg in np.linspace(0.0, 50.0, 26):  # step 2 degrees
        theta = math.radians(theta_deg)
        for M in np.linspace(-0.05, 0.05, 21):   # step 0.005
            for X in np.linspace(0.0, 100.0, 26):  # step 4
                L = loss(points, theta, M, X)
                if L < best[0]:
                    best = (L, theta, M, X)
    return best  # (L, theta, M, X)

# ----------------------------
# Local coordinate refinement
# (no SciPy; simple shrinking step search)
# ----------------------------
def refine(points: np.ndarray, theta: float, M: float, X: float,
           iters: int = 8,
           step_theta_deg: float = 2.0,
           step_M: float = 0.005,
           step_X: float = 4.0) -> Tuple[float, float, float, float]:

    Lbest = loss(points, theta, M, X)

    for _ in range(iters):
        improved = False

        # search theta
        for d in (-step_theta_deg, 0.0, step_theta_deg):
            th = theta + math.radians(d)
            th = max(0.0, min(math.radians(50.0), th))
            L = loss(points, th, M, X)
            if L < Lbest:
                theta, Lbest = th, L
                improved = True

        # search M
        for d in (-step_M, 0.0, step_M):
            m = M + d
            m = max(-0.05, min(0.05, m))
            L = loss(points, theta, m, X)
            if L < Lbest:
                M, Lbest = m, L
                improved = True

        # search X
        for d in (-step_X, 0.0, step_X):
            x = X + d
            x = max(0.0, min(100.0, x))
            L = loss(points, theta, M, x)
            if L < Lbest:
                X, Lbest = x, L
                improved = True

        # shrink steps each round
        step_theta_deg *= 0.5
        step_M *= 0.5
        step_X *= 0.5

        if not improved:
            # no improvement at this resolution; continue shrinking anyway
            pass

    return Lbest, theta, M, X

# ----------------------------
# Main
# ----------------------------
def main():
    if len(sys.argv) < 2:
        print("Usage: python fit_curve.py xy_data.csv")
        sys.exit(1)

    csv_path = sys.argv[1]
    points = load_xy(csv_path)

    print(f"Loaded {points.shape[0]} points.")

    print("Coarse search...")
    L0, theta0, M0, X0 = coarse_search(points)
    print(f"Coarse best: loss={L0:.6f}, theta={math.degrees(theta0):.3f}°, M={M0:.6f}, X={X0:.3f}")

    print("Refining...")
    Lf, theta, M, X = refine(points, theta0, M0, X0)
    print("\nFinal parameters:")
    print(f"  theta = {theta:.8f} rad  ({math.degrees(theta):.5f} deg)")
    print(f"  M     = {M:.8f}")
    print(f"  X     = {X:.8f}")
    print(f"  final loss (median L1 + penalty) = {Lf:.8f}")

    # Build the Desmos / latex-style string
    desmos = (
        f"(t*cos({theta:.8f}) - e^({M:.8f}*abs(t))*sin(0.3*t)*sin({theta:.8f}) + {X:.8f}, "
        f"42 + t*sin({theta:.8f}) + e^({M:.8f}*abs(t))*sin(0.3*t)*cos({theta:.8f}))"
    )
    with open("desmos_string.txt", "w") as f:
        f.write(desmos + "\n")

    print("\nDesmos submission string written to desmos_string.txt")
    print(desmos)

    t_est, s_est = rotate_and_split(points, theta, X)
    print(f"t_est range: min={t_est.min():.3f}, max={t_est.max():.3f}")

    # plot raw points
    plt.scatter(points[:,0], points[:,1], s=6, alpha=0.5, label='data')
    # plot fitted curve over t in [6,60]
    ts = np.linspace(6,60,600)
    xs = X + ts*np.cos(theta) - np.exp(M*np.abs(ts))*np.sin(0.3*ts)*np.sin(theta)
    ys = 42 + ts*np.sin(theta) + np.exp(M*np.abs(ts))*np.sin(0.3*ts)*np.cos(theta)
    plt.plot(xs, ys, linewidth=2, label='fit')
    plt.legend(); plt.xlabel('x'); plt.ylabel('y'); plt.title('Data vs Fitted Curve')
    plt.tight_layout(); plt.savefig('fit.png', dpi=150)
    print("Saved plot to fit.png")

if __name__ == "__main__":
    main()
