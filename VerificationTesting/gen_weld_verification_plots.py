#!/usr/bin/env python3
"""Four-panel verification plot for the DECCO welded rigid-body collision.

This implementation matches the definitions in the manuscript appendix:

(a) Projection of each body's linear momentum and the total momentum onto
    the INITIAL relative center-of-mass velocity direction.
(b) epsilon_P(t) = ||sum_i m_i v_i(t) - sum_i m_i v_i(0)||_2
                    / sum_i m_i |v_i(0)|.
(c) epsilon_L(t) = ||sum_i L_i(t) - sum_i L_i(0)||_2
                    / sum_i ||L_i(0)||_2,
    L_i(t) = x_i(t) x [m_i v_i(t)] + I_i^mono omega_i(t).
(d) For each original rigid body g, max_{i,k in G_g} |d_ik(t) - d_ik(0)|
    / max_{i,k in G_g} d_ik(0).

The time t=0 in these expressions is the first SELECTED saved snapshot
(--first-line).  The code does NOT recenter individual snapshots.

Assumptions:
  * utils.get_all_data(...) returns (pos, vel, w, radius, mass, moi, groups).
  * utils.get_energy(...) returns the snapshot times as its first result.
  * All arrays use CGS units and the same inertial reference frame.
  * Saved monomer velocities include rigid-body rotational contributions.
  * 'moi' is each monomer's intrinsic scalar moment of inertia.
  * Monomer order and original rigid-body group IDs remain fixed.
  * Use unwrapped coordinates if periodic boundaries are active.

Example:
    python validate_rigid_bodies.py PATH/TO/DATA \\
        --output verification_plot.png --save-data verification.csv --no-show

The default rigidity calculation includes ALL monomer pairs. Optional
--pair-limit N samples at most N pairs per body and is therefore approximate;
leave this option unset for the exact appendix definition.
"""

import argparse
import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.ticker import FormatStrFormatter

# This file is normally placed in DECCO/VerificationTesting/.
project_path = Path(__file__).resolve().parent.parent
sys.path.append(str(project_path / "utilities"))
import utils as u


def as_vectors(value, n=None, name="array"):
    arr = np.asarray(value, dtype=float)
    if arr.ndim != 2:
        raise ValueError(f"{name} must be 2D; received shape {arr.shape}")
    if arr.shape[1] == 3 and (n is None or arr.shape[0] == n):
        return arr
    if arr.shape[0] == 3 and (n is None or arr.shape[1] == n):
        return arr.T
    raise ValueError(f"{name} must have shape (N, 3) or (3, N); got {arr.shape}")


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
        raise ValueError("group IDs must be numeric")
    if not np.all(np.isfinite(ids)) or not np.all(ids == np.floor(ids)):
        raise ValueError("group IDs must be finite integers")
    return ids.astype(np.int64)


def get_snapshot(folder, data_index, line, relax):
    pos, vel, w, _radius, mass, moi, groups = u.get_all_data(
        folder, data_index=data_index, linenum=line, relax=relax
    )
    pos = as_vectors(pos, name="pos")
    n = len(pos)
    return (
        pos,
        as_vectors(vel, n=n, name="vel"),
        as_vectors(w, n=n, name="w"),
        as_scalars(mass, n, "mass"),
        as_scalars(moi, n, "moi"),
        as_group_ids(groups, n),
    )


def individual_angular_momenta(pos, vel, w, mass, moi):
    """Orbital + intrinsic spin angular momentum of every monomer (N, 3)."""
    return np.cross(pos, mass[:, None] * vel) + moi[:, None] * w


def pairs_for_bodies(group_ids, pair_limit, seed):
    """Return ORIGINAL monomer-index pairs, keyed by stable rigid-body ID."""
    rng = np.random.default_rng(seed)
    result = {}
    for body_id in np.unique(group_ids[group_ids >= 0]):
        members = np.flatnonzero(group_ids == body_id)
        ii, jj = np.triu_indices(len(members), k=1)
        if pair_limit and len(ii) > pair_limit:
            chosen = rng.choice(len(ii), size=pair_limit, replace=False)
            ii, jj = ii[chosen], jj[chosen]
        result[int(body_id)] = (members[ii], members[jj])
    return result


def main():
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    parser.add_argument("data_folder", help="Folder passed to utils.get_all_data")
    parser.add_argument("--steps", type=int, default=-1,
                        help="Number of saved snapshots (-1: all available)")
    parser.add_argument("--first-line", type=int, default=0)
    parser.add_argument("--data-index", type=int, default=-1)
    parser.add_argument("--relax", action="store_true")
    parser.add_argument("--pair-limit", type=int, default=0,
                        help="0: exact (all pairs, default); N>0: approximate sampled rigidity")
    parser.add_argument("--seed", type=int, default=1234,
                        help="Random seed, only used when sampling monomer pairs")
    parser.add_argument("--output", default="verification_plot.png")
    parser.add_argument("--save-data", help="Optional numerical CSV output path")
    parser.add_argument("--no-show", action="store_true")
    args = parser.parse_args()

    if args.first_line < 0 or args.pair_limit < 0:
        parser.error("--first-line and --pair-limit must be nonnegative")
    if args.steps != -1 and args.steps < 2:
        parser.error("--steps must be >= 2, or -1 for all saved snapshots")

    folder = str(Path(args.data_folder)) + "/"
    energy_time, *_ = u.get_energy(folder, data_index=args.data_index,
                                    relax=args.relax)
    energy_time = np.asarray(energy_time, dtype=float).reshape(-1)
    if not np.all(np.isfinite(energy_time)) or np.any(np.diff(energy_time) <= 0):
        raise ValueError("Energy times must be finite and strictly increasing")
    available = len(energy_time) - args.first_line
    if args.steps == -1:
        args.steps = available
    if not 2 <= args.steps <= available:
        raise ValueError(f"Requested {args.steps} snapshots starting at line "
                         f"{args.first_line}, with {len(energy_time)} timestamps")

    # Use the first selected snapshot as the t=0 baseline for all diagnostics.
    pos0, vel0, w0, mass0, moi0, groups0 = get_snapshot(
        folder, args.data_index, args.first_line, args.relax
    )
    body_ids = sorted(int(bid) for bid in np.unique(groups0) if bid >= 0)
    if len(body_ids) != 2:
        raise ValueError("Panel (a) defines the initial collision axis using two "
                         f"rigid bodies, but found group IDs {body_ids}")
    masks = {bid: groups0 == bid for bid in body_ids}

    # Project onto the relative, INITIAL center-of-mass velocity direction.
    # With body_ids [0, 1], this points from V_COM,1 toward V_COM,0, so the
    # incoming momentum of body 0 is positive in the supplied example plot.
    com_velocities = []
    for bid in body_ids:
        mask = masks[bid]
        com_velocities.append(
            np.sum(mass0[mask, None] * vel0[mask], axis=0) / np.sum(mass0[mask])
        )
    v_relative0 = com_velocities[0] - com_velocities[1]
    v_relative_norm = np.linalg.norm(v_relative0)
    if not np.isfinite(v_relative_norm) or v_relative_norm == 0:
        raise ValueError("Initial relative center-of-mass velocity is zero; "
                         "cannot define the collision direction")
    collision_axis = v_relative0 / v_relative_norm

    # Exactly the denominators written in the appendix (not ||P_total(0)||
    # or ||L_total(0)||, both of which may vanish in a COM-centered run).
    p_scale = np.sum(mass0 * np.linalg.norm(vel0, axis=1))
    l_i0 = individual_angular_momenta(pos0, vel0, w0, mass0, moi0)
    l_scale = np.sum(np.linalg.norm(l_i0, axis=1))
    if p_scale <= 0 or l_scale <= 0:
        raise ValueError("At least one appendix normalization is zero: "
                         f"p_scale={p_scale}, l_scale={l_scale}")

    pairs = pairs_for_bodies(groups0, args.pair_limit, args.seed)
    if args.pair_limit:
        print("Note: --pair-limit makes panel (d) an approximation to the "
              "all-pairs expression in the appendix.")
    d0 = {}
    dmax0 = {}
    for bid, (ii, jj) in pairs.items():
        if len(ii) == 0:
            raise ValueError(f"Rigid body {bid} has fewer than two monomers")
        d0[bid] = np.linalg.norm(pos0[ii] - pos0[jj], axis=1)
        dmax0[bid] = np.max(d0[bid])
        if dmax0[bid] == 0:
            raise ValueError(f"Rigid body {bid} has zero initial spatial extent")

    total_p = []
    total_l = []
    projected_p = {bid: [] for bid in body_ids}
    rigidity = {bid: [] for bid in body_ids}

    for k in range(args.steps):
        line = args.first_line + k
        try:
            pos, vel, w, mass, moi, groups = get_snapshot(
                folder, args.data_index, line, args.relax
            )
        except Exception as exc:
            raise RuntimeError(f"Cannot read snapshot line {line} from {folder}") from exc
        if len(pos) != len(pos0) or not np.array_equal(groups, groups0):
            raise ValueError(f"Monomer ordering/group membership changed at line {line}")
        if not np.array_equal(mass, mass0) or not np.array_equal(moi, moi0):
            raise ValueError(f"Monomer masses or intrinsic moments changed at line {line}")

        p_i = mass[:, None] * vel
        total_p.append(p_i.sum(axis=0))
        total_l.append(individual_angular_momenta(pos, vel, w, mass, moi).sum(axis=0))
        for bid in body_ids:
            projected_p[bid].append(p_i[masks[bid]].sum(axis=0) @ collision_axis)
            ii, jj = pairs[bid]
            d_t = np.linalg.norm(pos[ii] - pos[jj], axis=1)
            rigidity[bid].append(np.max(np.abs(d_t - d0[bid])) / dmax0[bid])

    p_arr = np.asarray(total_p)
    l_arr = np.asarray(total_l)
    dp = p_arr - p_arr[0]
    dl = l_arr - l_arr[0]
    p_error = np.linalg.norm(dp, axis=1) / p_scale
    l_error = np.linalg.norm(dl, axis=1) / l_scale
    total_p_parallel = p_arr @ collision_axis
    x = energy_time[args.first_line:args.first_line + args.steps]

    plt.rcParams.update({
        "font.size": 16,
        # "text.usetex": usetex,
        # "text.latex.preamble": r"\usepackage{amsmath} \usepackage{bm}" if usetex else "",
    })

    fig, axes = plt.subplots(
        4, 1, figsize=(8.5, 11), sharex=True, constrained_layout=True
    )
    for bid in body_ids:
        axes[0].plot(x, projected_p[bid], label=f"Body {bid}", linewidth=1.15)
    axes[0].plot(x, total_p_parallel, "k--", label="Total system", linewidth=1.4)
    axes[0].set_ylabel(r"$P_{\parallel}$ (g cm s$^{-1}$)")
    axes[0].legend(loc="lower right", ncol=3, fontsize=16)

    axes[1].plot(x, p_error, linewidth=1.15)
    axes[1].set_ylabel(r"$\epsilon_P(t)$")

    axes[2].plot(x, l_error, linewidth=1.15)
    axes[2].set_ylabel(r"$\epsilon_L(t)$")

    for bid in body_ids:
        axes[3].plot(x, rigidity[bid], label=f"Body {bid}", linewidth=1.1)
    axes[3].set_ylabel(r"$\Delta d_{\mathrm{max}}(t)/d_{\mathrm{max}}(0)$")
    axes[3].set_xlabel("Time (s)")
    axes[3].legend(loc="lower right", ncol=2, fontsize=16)

    for letter, ax in zip("abcd", axes):
        ax.text(0.01, 0.96, f"({letter})", transform=ax.transAxes,
                ha="left", va="top", fontsize=11)
        ax.grid(alpha=0.22)
        ax.ticklabel_format(axis="y", style="sci", scilimits=(-2, 2),
                            useMathText=True)
    axes[3].xaxis.set_major_formatter(FormatStrFormatter("%.5f"))

    out = Path(args.output)
    out.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out, dpi=200)

    print(f"Plot saved to {out}")
    print("Initial collision direction:", collision_axis)
    print(f"Linear momentum normalization:  {p_scale:.8e} g cm s^-1")
    print(f"Angular momentum normalization: {l_scale:.8e} g cm^2 s^-1")
    print(f"Max epsilon_P(t): {np.max(p_error):.6e}")
    print(f"Max epsilon_L(t): {np.max(l_error):.6e}")
    for bid in body_ids:
        print(f"Body {bid}: max normalized separation drift: "
              f"{max(rigidity[bid]):.6e}")

    if args.save_data:
        columns = [x, *dp.T, *dl.T, p_error, l_error]
        headers = ["time_s", "dPx", "dPy", "dPz", "dLx", "dLy", "dLz",
                   "epsilon_P", "epsilon_L"]
        for bid in body_ids:
            columns.append(np.asarray(projected_p[bid]))
            headers.append(f"P_parallel_body_{bid}")
        columns.append(total_p_parallel)
        headers.append("P_parallel_total")
        for bid in body_ids:
            columns.append(np.asarray(rigidity[bid]))
            headers.append(f"rigidity_body_{bid}")
        csv = Path(args.save_data)
        csv.parent.mkdir(parents=True, exist_ok=True)
        np.savetxt(csv, np.column_stack(columns), delimiter=",",
                   header=",".join(headers), comments="")
        print(f"Numerical data saved to {csv}")

    if not args.no_show:
        plt.show()
    plt.close(fig)


if __name__ == "__main__":
    main()