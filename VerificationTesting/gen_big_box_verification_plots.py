
#!/usr/bin/env python3
"""Validate a periodic DECCO Big Box coagulation run.

Uses the same utility reader as validate_rigid_bodies.py:
    u.get_all_data(data_folder, data_index=-1, linenum=-1, relax=False)
returning (pos, vel, w, radius, mass, moi, groups).

Example (plots all available snapshots by default):
    python validate_big_box.py PATH/TO/DATA \
        --save-data big_box_validation.csv --output big_box_validation.png

Requires numpy, matplotlib, and scipy (for efficient periodic neighbor queries).

Assumptions / limitations:
  * Energy-file timestamps from u.get_energy() correspond one-to-one with
    saved particle snapshots in u.get_all_data().
  * Box lengths are read from default_input.json via u.get_box_size()
    (comma-separated string, e.g. "0.01,0.01,0.01").
  * Positions are WRAPPED into [-box-lengths/2, +box-lengths/2) in
    each periodic dimension (box centered at the origin).
  * Monomers retain the same array ordering and identity across snapshots.
  * A displacement of > half the box length between two saved frames is
    interpreted as a boundary crossing. This inference is only reliable when
    each monomer travels < half a box length per axis between saved frames.
    Consecutive snapshots alone cannot prove how many times a particle wrapped.
  * "Periodic overlaps" are geometric overlaps across minimum-image boundaries
    at saved times, NOT direct evidence that DECCO actually applied contact forces.
    Collisions that occur between saved frames may be missed.
  * "groups" are aggregate IDs; negative IDs mean singleton monomers, each
    counted separately. Aggregate IDs may change as coagulation proceeds.
  * The full system is assumed to experience no external forces. Angular
    momentum is not used as a periodic-cell conservation diagnostic.
  * If particles can be reordered or removed/added, the identity-based checks
    and inferred boundary-crossing counters must be adapted first.
"""

import argparse
import os
from pathlib import Path
import sys

import matplotlib.pyplot as plt
import numpy as np
from scipy.spatial import cKDTree

from pathlib import Path

project_path = Path(__file__).resolve().parent.parent
sys.path.append(str(project_path / "utilities"))

import utils as u


def as_vectors(value, n=None, name="array"):
    arr = np.asarray(value, dtype=float)
    if arr.ndim != 2:
        raise ValueError(f"{name} must be a 2D array; got {arr.shape}")
    if arr.shape[1] == 3 and (n is None or arr.shape[0] == n):
        result = arr
    elif arr.shape[0] == 3 and (n is None or arr.shape[1] == n):
        result = arr.T
    else:
        raise ValueError(f"{name} must be shaped (N,3) or (3,N); got {arr.shape}")
    if not np.all(np.isfinite(result)):
        raise ValueError(f"{name} has non-finite values")
    return result


def as_scalars(value, n, name):
    arr = np.asarray(value, dtype=float).reshape(-1)
    if arr.size != n or not np.all(np.isfinite(arr)):
        raise ValueError(f"{name}: expected {n} finite values, got shape {arr.shape}")
    return arr


def as_group_ids(value, n):
    ids = np.asarray(value).reshape(-1)
    if ids.size != n or not np.issubdtype(ids.dtype, np.number):
        raise ValueError(f"groups: expected {n} numerical IDs")
    if not np.all(np.isfinite(ids)) or not np.all(ids == np.floor(ids)):
        raise ValueError("groups must contain finite integer IDs")
    return ids.astype(np.int64)


def aggregate_count(groups):
    # Each negative-ID monomer is treated as its own singleton, even if
    # multiple monomers share the same negative sentinel (e.g., -1).
    return int(np.unique(groups[groups >= 0]).size + np.count_nonzero(groups < 0))


def overlap_counts(pos, radii, lengths, origin, tolerance):
    """Count minimum-image overlaps, and the subset crossing periodic faces.

    Contacts are candidate geometric overlaps only. cKDTree searches within
    twice the largest radius, then tests each pair's actual sum of radii.
    """
    coords = np.mod(pos - origin, lengths)  # cKDTree periodic coords in [0,L)
    tree = cKDTree(coords, boxsize=lengths)
    search_radius = 2.0 * float(np.max(radii)) + tolerance
    try:
        pairs = tree.query_pairs(search_radius, output_type="ndarray")
    except TypeError:  # older SciPy
        pairs = np.asarray(list(tree.query_pairs(search_radius)), dtype=int).reshape(-1, 2)
    if pairs.size == 0:
        return 0, 0
    ii, jj = pairs[:, 0], pairs[:, 1]
    raw = coords[ii] - coords[jj]
    minimum = raw - lengths * np.rint(raw / lengths)
    overlap = np.einsum("ij,ij->i", minimum, minimum) <= (radii[ii] + radii[jj] + tolerance) ** 2
    crosses_face = np.any(np.abs(raw) > lengths / 2.0, axis=1)
    return int(np.count_nonzero(overlap)), int(np.count_nonzero(overlap & crosses_face))


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("data_folder", help="Folder passed directly to u.get_all_data")
    parser.add_argument("--steps", type=int, default=-1,
                        help="Number of saved snapshots to read (-1: all energy timestamps from --first-line)")
    parser.add_argument("--first-line", type=int, default=0)
    parser.add_argument("--data-index", type=int, default=-1)
    parser.add_argument("--relax", action="store_true")
    parser.add_argument("--contacts-every", type=int, default=1,
                        help="Evaluate minimum-image overlaps every this many saved snapshots")
    parser.add_argument("--contact-tolerance", type=float, default=0.0,
                        help="Nonnegative absolute position tolerance for geometric overlaps")
    parser.add_argument("--output", default="big_box_validation.png")
    parser.add_argument("--save-data", help="Optional numeric CSV path")
    parser.add_argument("--no-show", action="store_true", help="Do not display the plot")
    args = parser.parse_args()
    if args.first_line < 0 or (args.steps != -1 and args.steps < 2) or args.contacts_every < 1:
        parser.error("--first-line must be >= 0; --steps must be >= 2 or -1; --contacts-every must be >= 1")
    if args.contact_tolerance < 0:
        parser.error("--contact-tolerance must be nonnegative")
    folder = args.data_folder.rstrip("/") + "/"
    # Read periodic box dimensions from default_input.json in this run folder.
    box_dims = u.get_box_size(folder)
    # default_input.json stores boxDims as e.g. "0.01,0.01,0.01".
    if isinstance(box_dims, str):
        box_dims = [part.strip() for part in box_dims.split(",")]
    try:
        lengths = np.asarray(box_dims, dtype=float)
    except (ValueError, TypeError) as exc:
        parser.error(f"Could not parse boxDims from input.json: {exc}")
    if lengths.shape != (3,) or not np.all(np.isfinite(lengths)) or np.any(lengths <= 0):
        parser.error("default_input.json boxDims must contain three finite positive box lengths")
    origin = -0.5 * lengths  # Periodic box centered at (0, 0, 0)
    print(f"Periodic box dimensions from {folder}input.json: {lengths}")
    # DECCO writes the physical snapshot times to its energy file.
    # As in gen_verification_plots_time, require matching row order with simData.
    energy_time, _, _, _, _, _ = u.get_energy(
        folder, data_index=args.data_index, relax=args.relax
    )
    energy_time = np.asarray(energy_time, dtype=float).reshape(-1)
    if not np.all(np.isfinite(energy_time)):
        raise ValueError("Energy file contains non-finite timestamps")
    if np.any(np.diff(energy_time) <= 0):
        raise ValueError("Energy timestamps must be strictly increasing")
    available = len(energy_time) - args.first_line
    if args.steps == -1:
        args.steps = available
    if args.steps < 1 or args.steps > available:
        raise ValueError(
            f"Requested {args.steps} snapshots starting at line {args.first_line}, "
            f"but energy data has {len(energy_time)} rows. "
            "Check that energy rows and particle snapshots correspond one-to-one."
        )
    x = energy_time[args.first_line:args.first_line + args.steps]
    xlabel = "Time (s)"
    p_series, n_series, m_series, aggregate_series = [], [], [], []
    crossed_series, overlap_series, periodic_overlap_series = [], [], []
    previous_pos = None
    initial_mass = None
    initial_radius = None
    initial_count = None
    cumulative_crossings = 0
    maximum_frame_displacement_ratio = 0.0

    for k in range(args.steps):
        line = args.first_line + k
        try:
            pos_raw, vel_raw, _, radius_raw, mass_raw, _, groups_raw = u.get_all_data(
                folder, data_index=args.data_index, linenum=line, relax=args.relax
            )
        except Exception as exc:
            raise RuntimeError(f"Could not read saved line {line} from {folder}") from exc

        pos = as_vectors(pos_raw, name="pos")
        n = len(pos)
        if n == 0:
            raise ValueError(f"No monomers at saved line {line}")
        vel = as_vectors(vel_raw, n=n, name="vel")
        radii = as_scalars(radius_raw, n, "radius")
        mass = as_scalars(mass_raw, n, "mass")
        groups = as_group_ids(groups_raw, n)
        if np.any(radii <= 0) or np.any(mass <= 0):
            raise ValueError(f"Nonpositive monomer radius or mass at saved line {line}")

        if k == 0:
            initial_count = n
            initial_mass = mass.copy()
            initial_radius = radii.copy()
        elif n != initial_count:
            raise ValueError(f"Monomer count changed at saved line {line}: {initial_count} -> {n}")
        elif not np.array_equal(mass, initial_mass) or not np.array_equal(radii, initial_radius):
            raise ValueError(f"Monomer masses/radii or ordering changed at saved line {line}")

        # Allow tiny roundoff in the saved coordinates, but flag grossly
        # out-of-box positions rather than silently counting invalid crossings.
        box_slack = 1e-10 * lengths
        outside = np.any(
            (pos < origin - box_slack) |
            (pos > origin + lengths + box_slack),
            axis=1
        )

        bad_indices = np.where(outside)[0]

        if len(bad_indices) > 0:
            raise ValueError(
                f"Positions outside specified periodic box at saved line {line}.\n"
                f"Particle indices: {bad_indices}\n"
                f"Particle positions:\n{pos[bad_indices]}"
            )

        p_series.append(np.sum(mass[:, None] * vel, axis=0))
        n_series.append(n)
        m_series.append(float(np.sum(mass)))
        aggregate_series.append(aggregate_count(groups))

        if previous_pos is not None:
            displacement = pos - previous_pos
            # Detect discontinuous coordinate jumps, not actual collision events.
            crossings = int(np.count_nonzero(np.abs(displacement) > lengths / 2))
            cumulative_crossings += crossings
            maximum_frame_displacement_ratio = max(
                maximum_frame_displacement_ratio, float(np.max(np.abs(displacement) / lengths))
            )
        crossed_series.append(cumulative_crossings)
        previous_pos = pos.copy()

        if k % args.contacts_every == 0:
            total_overlaps, periodic_overlaps = overlap_counts(
                pos, radii, lengths, origin, args.contact_tolerance
            )
        else:
            total_overlaps, periodic_overlaps = np.nan, np.nan
        overlap_series.append(total_overlaps)
        periodic_overlap_series.append(periodic_overlaps)

    momenta = np.asarray(p_series)
    counts = np.asarray(n_series)
    masses = np.asarray(m_series)
    aggregates = np.asarray(aggregate_series)
    crossings = np.asarray(crossed_series)
    overlaps = np.asarray(overlap_series, dtype=float)
    periodic_overlaps = np.asarray(periodic_overlap_series, dtype=float)

    # This scale is meaningful even if the net initial momentum is near zero.
    p0_raw, v0_raw, _, _, m0_raw, _, _ = u.get_all_data(
        folder, data_index=args.data_index, linenum=args.first_line, relax=args.relax
    )
    v0 = as_vectors(v0_raw, n=initial_count, name="initial velocities")
    m0 = as_scalars(m0_raw, initial_count, "initial masses")
    p_scale = max(np.linalg.norm(momenta[0]), np.sum(m0 * np.linalg.norm(v0, axis=1)),
                  np.finfo(float).tiny)
    p_delta = momenta - momenta[0]
    relative_p = np.linalg.norm(p_delta, axis=1) / p_scale
    relative_mass = np.abs(masses - masses[0]) / masses[0]

    fig, axes = plt.subplots(3, 2, figsize=(12, 10), sharex=True, constrained_layout=True)
    ax = axes[0, 0]
    for j, component in enumerate("xyz"):
        ax.plot(x, momenta[:, j], label=f"$P_{component}$", linewidth=1.4)
    ax.set_ylabel("Total linear momentum")
    ax.legend()

    axes[0, 1].plot(x, relative_p, linewidth=1.6)
    axes[0, 1].set_ylabel("Relative linear momentum error")

    axes[1, 0].plot(x, counts, label="Monomers", linewidth=1.6)
    axes[1, 0].plot(x, aggregates, label="Aggregates", linewidth=1.6)
    axes[1, 0].set_ylabel("Number")
    axes[1, 0].legend()

    axes[1, 1].plot(x, relative_mass, linewidth=1.6)
    axes[1, 1].set_ylabel("Relative total mass error")

    axes[2, 0].step(x, crossings, where="post", linewidth=1.6)
    axes[2, 0].set_ylabel("Inferred cumulative boundary crossings")

    axes[2, 1].plot(x, overlaps, label="All overlaps", linewidth=1.4)
    axes[2, 1].plot(x, periodic_overlaps, label="Across boundary", linewidth=1.4)
    axes[2, 1].set_ylabel("Geometric overlapping pairs")
    axes[2, 1].legend()

    for ax in axes.flat:
        ax.grid(alpha=0.3)
        ax.ticklabel_format(axis="y", style="sci", scilimits=(-3, 3))
    for ax in axes[-1, :]:
        ax.set_xlabel(xlabel)
    fig.suptitle("Periodic Big Box: conservation and boundary diagnostics")

    out = Path(args.output)
    out.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out, dpi=200)
    print(f"Plot saved to {out}")
    print(f"Max relative linear momentum change: {np.max(relative_p):.6e}")
    print(f"Max relative total mass change: {np.max(relative_mass):.6e}")
    print(f"Monomer count: {initial_count} (unchanged)")
    print(f"Aggregate count: {aggregates[0]} -> {aggregates[-1]}")
    print(f"Inferred boundary-crossing events: {crossings[-1]}")
    print(f"Largest sampled total overlap count: {int(np.nanmax(overlaps))}")
    print(f"Largest sampled cross-boundary overlap count: {int(np.nanmax(periodic_overlaps))}")
    if crossings[-1] == 0:
        print("WARNING: No wrapped-position jumps detected in saved snapshots.")
    if np.nanmax(periodic_overlaps) == 0:
        print("WARNING: No cross-boundary geometric overlaps found in sampled snapshots.")
    if maximum_frame_displacement_ratio > 0.5:
        print("CAUTION: Large between-snapshot coordinate jumps found. Boundary-crossing "
              "counts can be ambiguous, particularly if a monomer wraps multiple times.")
    if np.any(np.diff(aggregates) > 0):
        print("NOTE: Aggregate count increased in at least one interval (could reflect fragmentation "
              "or bookkeeping behavior).")

    if args.save_data:
        csv = Path(args.save_data)
        csv.parent.mkdir(parents=True, exist_ok=True)
        data = np.column_stack([
            x, *momenta.T, *p_delta.T, relative_p, counts, masses,
            relative_mass, aggregates, crossings, overlaps, periodic_overlaps
        ])
        headers = ["time_s", "Px", "Py", "Pz", "dPx", "dPy", "dPz",
                   "relative_dP", "monomer_count", "total_mass", "relative_mass_error",
                   "aggregate_count", "inferred_cumulative_boundary_crossings",
                   "geometric_overlaps", "geometric_cross_boundary_overlaps"]
        np.savetxt(csv, data, delimiter=",", header=",".join(headers), comments="")
        print(f"Numeric data saved to {csv}")
    if not args.no_show:
        plt.show()
    plt.close(fig)


if __name__ == "__main__":
    main()
