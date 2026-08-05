"""
Generate the three Parker-spiral input files the Strauss legacy `code.f90`
expects:

    Parker_CosPsi_108_Phi.dat   cos(ψ)  on (N=121, M=108) grid, row-major i then j
    Parker_SinPsi_108_Phi.dat   sin(ψ)  same grid
    Parker_Focus_108.dat        1/L     same grid (the legacy code stores L = 1/placeholder)

Grid (matches the hard-coded constants in code.f90):
    r_i = 0.05, ... 1.20 AU,  N = 121, Δr = (1.20-0.05)/(N-1) = 0.0095833 AU
    φ_j = -π, ... π-Δφ,        M = 108,  Δφ = 2π/M

Parker spiral with constant solar wind speed V_sw:
    tan(ψ) = Ω r / V_sw          (with helio-equatorial sin(θ)=1)
    cos(ψ) = 1 / √(1+tan²ψ)
    sin(ψ) = tan(ψ) · cos(ψ)
    1/L    = focusing length inverse;  for an unmagnetised B = B_0 (r0/r)² along
             the spiral, L⁻¹ = -∂_s ln B² / 2  (s along the spiral coordinate)
    Equivalently for the Parker spiral with strict r-φ geometry:
            1/L = (cos²ψ / r) · [2 - (sin²ψ - cos²ψ)/cos²ψ · (Ω r / V_sw)²·... ]
    For brevity we use the standard Parker form (Roelof 1969, Kocharov 1998):
        1/L = (1/B) dB/ds   with   B(r) = B_0 (r₀/r)² √(1 + (Ω r/V_sw)²)
        ds/dr = 1/cos(ψ)
        => 1/L = cos(ψ) · d/dr [ln B(r)]
                = cos(ψ) · [-2/r - (Ω/V_sw)² r / (1 + (Ω r/V_sw)²)]
        => focusing length L > 0 because dlnB/ds < 0 outward; we store |1/L|.

Usage:
    python make_parker_spiral.py --vsw 400 --outdir ~/work/sep2d/legacy

The output files are plain ASCII, one value per line, in legacy Fortran read
order: outer loop over i (radial), inner over j (azimuthal).
"""

from __future__ import annotations

import argparse
import os
from pathlib import Path

import numpy as np


# Grid constants — default matches code.f90's native PARAMETER grid.
# Overridable via --N/--M so the benchmark harness can generate Parker
# inputs for a recompiled-at-different-resolution legacy binary (the
# cpu_vs_gpu scaling curve). NK (pitch angle) does NOT appear here —
# cos ψ / sin ψ / 1/L are position-only, so NK can vary freely.
N = 121
M = 108
R_MIN = 0.05      # AU
R_MAX = 1.20      # AU

# Solar rotation rate (sidereal): Ω = 2π / (25.38 d)   ≈ 2.866e-6 rad/s
OMEGA_RAD_PER_S = 2.0 * np.pi / (25.38 * 86400.0)

# Convert Ω to AU·hr units used by code.f90 (speed_of_light = 7.2 AU/hr).
# 1 AU = 1.495978707e11 m; 1 hr = 3600 s
AU_M = 1.495978707e11
HR_S = 3600.0

# We need (Ω r / V_sw) dimensionless, so do everything in SI then ratio.


def parker(vsw_kms: float, n: int = N, m: int = M, focusing: str = "legacy",
           r_min: float = R_MIN, r_max: float = R_MAX
           ) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Parker-spiral input fields for the legacy Fortran code.

    ``focusing`` selects the sign of the spiral term in 1/L:
      * ``"legacy"``    — 1/L = cosψ·[2/r **+** (Ω/V)²r/(1+tan²ψ)], the form
        this generator has always emitted and that the published legacy runs
        used. It OVER-focuses (L ≈ 1.7× too short at 1 AU).
      * ``"corrected"`` — 1/L = cosψ·[2/r **−** (Ω/V)²r/(1+tan²ψ)], matching
        sep2d_hpc/parker.py (fixed 2026-06-03). The spiral makes B fall off
        more slowly outward, which WEAKENS focusing; the large-r limit B∝1/r
        gives 1/L→cosψ/r and textbook L(1 AU)≈1 AU, both only reproduced by
        the minus form.

    Because 1/L reaches the legacy code as an input file, this switch lets the
    legacy binary be run with either focusing without touching the Fortran —
    which is what makes an honest legacy↔JAX comparison possible.
    """
    if focusing not in ("legacy", "corrected"):
        raise ValueError(f"focusing must be legacy|corrected, got {focusing!r}")
    vsw = vsw_kms * 1e3  # m/s
    r_au = np.linspace(r_min, r_max, n)            # (n,)
    r_m = r_au * AU_M
    # tan ψ = Ω r / V_sw (helio-equatorial, sin θ = 1)
    tanpsi = OMEGA_RAD_PER_S * r_m / vsw           # (n,)
    cospsi = 1.0 / np.sqrt(1.0 + tanpsi**2)
    sinpsi = tanpsi * cospsi
    # 1/L (units: 1/AU). cospsi is dimensionless; the spiral term is 1/m here.
    spiral = (OMEGA_RAD_PER_S / vsw) ** 2 * r_m / (1.0 + tanpsi**2)
    sign = +1.0 if focusing == "legacy" else -1.0
    inv_L_per_m = cospsi * (2.0 / r_m + sign * spiral)
    inv_L_per_au = inv_L_per_m * AU_M
    # The legacy code stores `placeholder = 1/L` and then sets L = 1/placeholder.
    # We supply 1/L (positive: focusing length is positive outward).
    # Tile to (n, m) — Parker spiral is φ-independent in this analytic form.
    cospsi_2d = np.broadcast_to(cospsi[:, None], (n, m)).copy()
    sinpsi_2d = np.broadcast_to(sinpsi[:, None], (n, m)).copy()
    inv_L_2d = np.broadcast_to(inv_L_per_au[:, None], (n, m)).copy()
    return cospsi_2d, sinpsi_2d, inv_L_2d


def write_dat(path: Path, arr: np.ndarray) -> None:
    """Fortran reads i (outer, radial), j (inner, azimuthal). One value per line."""
    with open(path, "w") as f:
        for i in range(arr.shape[0]):
            for j in range(arr.shape[1]):
                f.write(f"{arr[i, j]:.10E}\n")


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--vsw", type=float, default=400.0,
                    help="Solar wind speed [km/s] (default 400)")
    ap.add_argument("--N", type=int, default=N, help=f"radial grid pts (default {N})")
    ap.add_argument("--M", type=int, default=M, help=f"azimuthal grid pts (default {M})")
    ap.add_argument("--focusing", choices=("legacy", "corrected"), default="legacy",
                    help="sign of the spiral term in 1/L: 'legacy' (+, over-focusing, "
                         "as published) or 'corrected' (-, matches sep2d parker.py)")
    ap.add_argument("--r-max", type=float, default=R_MAX,
                    help=f"outer radius [AU] (default {R_MAX}; must match the "
                         f"legacy build's delta_r constant)")
    ap.add_argument("--outdir", type=Path, default=Path.cwd(),
                    help="Output directory (default cwd)")
    args = ap.parse_args()

    args.outdir.mkdir(parents=True, exist_ok=True)
    cospsi, sinpsi, inv_L = parker(args.vsw, n=args.N, m=args.M,
                                   focusing=args.focusing, r_max=args.r_max)

    # Filenames keep the legacy "108" tag (code.f90 reads these exact names);
    # the actual array size follows --N/--M, and the legacy READ loop reads
    # N×M values per its own PARAMETER grid, so the two must be kept in sync
    # by the caller (benchmark/cpu_vs_gpu.py patches both together).
    write_dat(args.outdir / "Parker_CosPsi_108_Phi.dat", cospsi)
    write_dat(args.outdir / "Parker_SinPsi_108_Phi.dat", sinpsi)
    write_dat(args.outdir / "Parker_Focus_108.dat", inv_L)

    print(f"V_sw         : {args.vsw} km/s")
    print(f"r grid       : {args.N} pts, {R_MIN}–{args.r_max} AU  ({args.focusing} focusing)")
    print(f"phi grid     : {args.M} pts, -π…π-Δφ")
    print(f"wrote 3 files to {args.outdir}/")


if __name__ == "__main__":
    main()
