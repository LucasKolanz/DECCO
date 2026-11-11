
#!/usr/bin/env python3
"""Validate conservation laws and rigidity in a DECCO rigid-body collision.

Requires a module named Utils importable as `import Utils as u` (override with
--module); it must implement:
    u.get_all_data(data_folder, data_index=-1, linenum=-1, relax=False)
returning (pos, vel, w, radius, mass, moi, groups).

Run from an environment where your DECCO Python utilities are importable.

Examples:
    python validate_rigid_bodies.py PATH/TO/DATA --steps 150
    python validate_rigid_bodies.py PATH/TO/DATA --steps 150 \
        --output validation.png --save-data validation.csv --no-show

Assumptions:
  - Each --linenum from 0 to --steps-1 is a saved snapshot, at fixed ordering.
  - Monomer positions and velocities refer to the same inertial reference frame.
  - Saved velocities include the rigid-body rotational velocity omega x r_rel.
  - "moi" is each monomer's scalar intrinsic moment of inertia, and "w" is
    its angular velocity; spin angular momentum is sum_i I_i * w_i.
  - Unwrapped positions should be used if periodic boundaries are enabled.
  - Group IDs identify the same rigid bodies in every snapshot. Negative group IDs
    indicate ungrouped monomers and are excluded from the rigidity check.

All positions can be translated by a constant origin; momentum conservation
plots use changes in total angular momentum, which are origin-independent if
linear momentum is conserved. Do not recenter each snapshot independently.
"""

import argparse
import importlib
import sys
import os
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

from pathlib import Path

project_path = Path(__file__).resolve().parent.parent
sys.path.append(str(project_path / "utilities"))

import utils as u


def as_vectors(value, n=None, name="array"):
    arr = np.asarray(value, dtype=float)
    if arr.ndim != 2:
        raise ValueError(f"{name} must be a 2D array; got shape {arr.shape}")
    if arr.shape[1] == 3 and (n is None or arr.shape[0] == n):
        return arr
    if arr.shape[0] == 3 and (n is None or arr.shape[1] == n):
        return arr.T
    raise ValueError(f"{name} must be shaped (N,3) or (3,N); got {arr.shape}")


def as_scalars(value, n, name):
    arr = np.asarray(value, dtype=float).reshape(-1)
    if arr.size != n:
        raise ValueError(f"{name}: expected {n} values, got {arr.size}")
    return arr


def as_group_ids(value, n):
    ids = np.asarray(value).reshape(-1)
    if ids.size != n:
        raise ValueError(f"groups: expected {n} values, got {ids.size}")
    if not np.issubdtype(ids.dtype, np.number):
        raise ValueError("groups must contain numerical rigid-body indices")
    if not np.all(np.isfinite(ids)) or not np.all(ids == np.floor(ids)):
        raise ValueError("groups must contain finite integer indices")
    return ids.astype(np.int64)


def pairs_for_bodies(body_ids, pair_limit, rng):
    """All pairs in small bodies; reproducibly sample up to pair_limit in large ones."""
    result = {}
    for body_id in np.unique(body_ids[body_ids >= 0]):
        members = np.flatnonzero(body_ids == body_id)
        ii, jj = np.triu_indices(len(members), k=1)
        if len(ii) > pair_limit:
            take = rng.choice(len(ii), size=pair_limit, replace=False)
            ii, jj = ii[take], jj[take]
        result[str(body_id)] = (members[ii], members[jj])
    return result


def compute_momenta(pos, vel, w, mass, moi, include_spin):
    p = (mass[:, None] * vel).sum(axis=0)
    orbital_l = np.cross(pos, mass[:, None] * vel).sum(axis=0)
    spin_l = (moi[:, None] * w).sum(axis=0) if include_spin else np.zeros(3)
    return p, orbital_l + spin_l


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("data_folder", help="Folder passed directly to u.get_all_data")
    parser.add_argument("--steps", type=int, default=-1, help="Number of saved snapshots (-1 = all available)")
    parser.add_argument("--first-line", type=int, default=0)
    parser.add_argument("--data-index", type=int, default=-1)
    parser.add_argument("--relax", action="store_true", default=False)
    parser.add_argument("--pair-limit", type=int, default=44850, help="Max checked monomer pairs per body")
    parser.add_argument("--seed", type=int, default=1234)
    parser.add_argument("--omit-spin", action="store_true", help="Omit intrinsic monomer spin from total angular momentum")
    parser.add_argument("--output", default="rigid_body_validation.png")
    parser.add_argument("--save-data", help="Optional path for numeric validation CSV")
    parser.add_argument("--no-show", action="store_true", help="Save results without opening a plot window")
    args = parser.parse_args()
    if args.first_line < 0 or args.pair_limit < 1:
        parser.error("--first-line must be >= 0 and --pair-limit >= 1")
    if args.steps != -1 and args.steps < 2:
        parser.error("--steps must be >= 2, or -1 for all saved outputs")

    if not args.data_folder.endswith("/"):
        args.data_folder += "/"

    energy_time, _, _, _, _, _ = u.get_energy(
        args.data_folder, data_index=args.data_index, relax=args.relax
    )
    energy_time = np.asarray(energy_time, dtype=float).reshape(-1)
    if not np.all(np.isfinite(energy_time)) or np.any(np.diff(energy_time) <= 0):
        raise ValueError("Energy timestamps must be finite and strictly increasing")
    available = len(energy_time) - args.first_line
    if args.steps == -1:
        args.steps = available
    if args.steps < 2 or args.steps > available:
        raise ValueError(
            f"Requested {args.steps} snapshots starting at {args.first_line}, "
            f"but there are {len(energy_time)} energy timestamps."
        )

    p_list, l_list = [], []
    rigidity = {}
    initial_distances = {}
    n0 = None
    body_ids = None
    for k in range(args.steps):
        line = args.first_line + k
        try:
            pos_raw, vel_raw, w_raw, radius_raw, mass_raw, moi_raw, groups_raw = u.get_all_data(
                args.data_folder, data_index=args.data_index, linenum=line, relax=args.relax
            )
        except Exception as exc:
            raise RuntimeError(f"Could not read saved line {line} from {args.data_folder}") from exc
        pos = as_vectors(pos_raw, name="pos")
        n = len(pos)
        vel = as_vectors(vel_raw, n=n, name="vel")
        w = as_vectors(w_raw, n=n, name="w")
        mass = as_scalars(mass_raw, n, "mass")
        moi = as_scalars(moi_raw, n, "moi")
        current_ids = as_group_ids(groups_raw, n)
        if k == 0:
            n0 = n
            body_ids = current_ids.copy()
            pairs = pairs_for_bodies(body_ids, args.pair_limit, np.random.default_rng(args.seed))
            for label, (ii, jj) in pairs.items():
                initial_distances[label] = np.linalg.norm(pos[ii] - pos[jj], axis=1)
                rigidity[label] = []
        elif n != n0:
            raise ValueError(f"Monomer count changed at saved line {line}: {n0} -> {n}")
        elif not np.array_equal(current_ids, body_ids):
            raise ValueError(f"Group membership changed at saved line {line}; "
                             "rigidity pairs require fixed membership and monomer ordering")

        p, l = compute_momenta(pos, vel, w, mass, moi, include_spin=not args.omit_spin)
        p_list.append(p)
        l_list.append(l)
        for label, (ii, jj) in pairs.items():
            d0 = initial_distances[label]
            if not len(d0):
                rigidity[label].append(0.0)
                continue
            distances = np.linalg.norm(pos[ii] - pos[jj], axis=1)
            # Normalize to the body's initial largest tested internal separation.
            scale = max(d0.max(), np.finfo(float).tiny)
            rigidity[label].append(np.max(np.abs(distances - d0)) / scale)

    p_arr = np.asarray(p_list)
    l_arr = np.asarray(l_list)
    # Fixed, nonzero scales also work when P(0) or L(0) is (nearly) zero.
    # A characteristic angular-momentum scale based on individual initial contributions.
    # Re-read the first snapshot to avoid basing the normalization on the last one.
    pos0, vel0, w0, _, mass0, moi0, _ = u.get_all_data(
        args.data_folder, data_index=args.data_index, linenum=args.first_line, relax=args.relax
    )
    pos0 = as_vectors(pos0, n=n0, name="pos0")
    vel0 = as_vectors(vel0, n=n0, name="vel0")
    w0 = as_vectors(w0, n=n0, name="w0")
    mass0 = as_scalars(mass0, n0, "mass0")
    moi0 = as_scalars(moi0, n0, "moi0")
    p_scale = max(np.linalg.norm(p_arr[0]), np.sum(mass0 * np.linalg.norm(vel0, axis=1)), np.finfo(float).tiny)
    l_parts = np.cross(pos0, mass0[:, None] * vel0)
    if not args.omit_spin:
        l_parts += moi0[:, None] * w0
    l_scale = max(np.linalg.norm(l_arr[0]), np.sum(np.linalg.norm(l_parts, axis=1)), np.finfo(float).tiny)

    dp = p_arr - p_arr[0]
    dl = l_arr - l_arr[0]
    p_error = np.linalg.norm(dp, axis=1) / p_scale
    l_error = np.linalg.norm(dl, axis=1) / l_scale
    x = energy_time[args.first_line:args.first_line + args.steps]
    xlabel = "Time (s)"

    fig, axes = plt.subplots(3, 1, figsize=(9, 10), sharex=True, constrained_layout=True)
    axes[0].plot(x, p_error, label="|P(t) - P(0)| / P_scale", linewidth=1.7)
    axes[0].set_ylabel("Relative linear momentum error")
    axes[1].plot(x, l_error, label="|L(t) - L(0)| / L_scale", linewidth=1.7)
    axes[1].set_ylabel("Relative angular momentum error")
    for label, errors in rigidity.items():
        axes[2].plot(x, errors, label=f"Body {label}")
    axes[2].set_ylabel("Max normalized pair-distance drift")
    axes[2].set_xlabel(xlabel)
    for ax in axes:
        ax.grid(alpha=0.3)
        ax.ticklabel_format(axis="y", style="sci", scilimits=(-3, 3))
    if rigidity:
        axes[2].legend()
    else:
        print("Warning: no nonnegative group IDs found; rigid-body constraint panel is empty")
    fig.suptitle("Rigid-body collision: conservation and constraints")

    out = Path(args.output)
    out.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out, dpi=200)

    print(f"Plot saved to {out}")
    print(f"Max relative linear momentum change: {np.max(p_error):.6e}")
    print(f"Max relative angular momentum change: {np.max(l_error):.6e}")
    for label, errors in rigidity.items():
        print(f"Body {label}: max normalized pair-distance drift = {max(errors):.6e}")
    if args.save_data:
        data = np.column_stack([x, *dp.T, *dl.T, p_error, l_error, *[rigidity[label] for label in rigidity]])
        headers = ["time_or_line", "dPx", "dPy", "dPz", "dLx", "dLy", "dLz", "relative_dP", "relative_dL"]
        headers += [f"rigidity_{label}" for label in rigidity]
        csv = Path(args.save_data)
        csv.parent.mkdir(parents=True, exist_ok=True)
        np.savetxt(csv, data, delimiter=",", header=",".join(headers), comments="")
        print(f"Numeric data saved to {csv}")
    if not args.no_show:
        plt.show()
    plt.close(fig)


if __name__ == "__main__":
    main()
