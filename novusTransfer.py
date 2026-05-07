"""
Transfer completed SpaceLab/DECCO jobs from the Novus cluster to the local machine.

This version intentionally uses the system OpenSSH commands (`ssh` and `scp`)
instead of Paramiko for authentication. That means it uses the exact same
~/.ssh/config behavior as this command:

    ssh Novus

This is useful when Paramiko fails with errors like:

    paramiko.ssh_exception.SSHException: No existing session

Run a dry-run first:

    python3 novusTransferNew.py --dry-run --job-set BAPAWELD

Then actually copy:

    python3 novusTransferNew.py --job-set BAPAWELD
"""

from __future__ import annotations

import argparse
import json
import os
import subprocess
import sys
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable


# -----------------------------------------------------------------------------
# Project imports
# -----------------------------------------------------------------------------

PROJECT_PATH = Path(__file__).resolve().parent
sys.path.append(str(PROJECT_PATH / "utilities"))

import utils as u  # noqa: E402


# -----------------------------------------------------------------------------
# Configuration
# -----------------------------------------------------------------------------

REMOTE_BASE_FOLDER = Path("/home/kolanzl/novus/kolanzl/SpaceLab_data")
DEFAULT_INPUT = PROJECT_PATH / "default_files" / "default_input.json"


@dataclass(frozen=True)
class TransferConfig:
    host_alias: str = "Novus"
    job_set_names: tuple[str, ...] = ("BAPAWELD",)
    attempts: range = range(30)
    m_values: tuple[int, ...] = (3, 5, 10, 15, 20, 30, 50, 60, 100)
    n_default: int = 300
    temps: tuple[int, ...] = (1000,)
    cbapa_c: int = 30
    dry_run: bool = False
    verbose_ssh: bool = False


def local_job_parent(job_set_name: str) -> str:
    """Choose the local parent folder used by older SpaceLab transfer scripts."""
    if job_set_name == "lognorm":
        return "jobsCosine"
    if job_set_name == "const":
        return "jobsNovus"
    return "jobs"


def n_for_job(job_set_name: str, m: int, n_default: int, cbapa_c: int) -> int:
    """CBAPA jobs use N = C*M; other current jobs use the configured default N."""
    if job_set_name == "CBAPA":
        return cbapa_c * m
    return n_default


# -----------------------------------------------------------------------------
# OpenSSH helpers
# -----------------------------------------------------------------------------

def run_command(command: list[str], *, check: bool = False) -> subprocess.CompletedProcess:
    """Run a shell command safely without shell=True."""
    return subprocess.run(
        command,
        text=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        check=check,
    )


def ssh_command(config: TransferConfig, remote_command: str) -> list[str]:
    """Build an ssh command that uses the user's normal ~/.ssh/config."""
    command = ["ssh"]

    if config.verbose_ssh:
        command.append("-v")

    # BatchMode prevents the script from hanging on password prompts.
    # ConnectTimeout keeps failures quick if the cluster/login node is down.
    command.extend([
        "-o", "BatchMode=yes",
        "-o", "ConnectTimeout=20",
        config.host_alias,
        remote_command,
    ])

    return command


def remote_file_exists(config: TransferConfig, remote_path: Path) -> bool:
    """
    Check whether a remote file exists using the system ssh command.

    Return code meanings:
      0 -> file exists
      1 -> ssh worked, but test -f returned false
      other -> ssh/connection/config problem
    """
    remote_command = f"test -f {shlex_quote(str(remote_path))}"
    result = run_command(ssh_command(config, remote_command))

    if result.returncode == 0:
        return True

    if result.returncode == 1:
        return False

    print("SSH command failed while checking remote file.")
    print("Command:", " ".join(ssh_command(config, remote_command)))
    print("stderr:")
    print(result.stderr.strip())
    raise RuntimeError(f"ssh failed with exit code {result.returncode}")


def scp_get_folder(config: TransferConfig, remote_folder: Path, local_parent: Path) -> bool:
    """
    Copy a remote folder recursively into local_parent using system scp.

    Example:
        remote_folder = /remote/jobs/BAPAWELD_0/M_3/N_300/T_1000/
        local_parent  = /local/jobs/BAPAWELD_0/M_3/N_300/

    This creates/updates:
        /local/jobs/BAPAWELD_0/M_3/N_300/T_1000/
    """
    local_parent.mkdir(parents=True, exist_ok=True)

    remote_spec = f"{config.host_alias}:{str(remote_folder)}"

    command = [
        "scp",
        "-r",
        "-o", "BatchMode=yes",
        "-o", "ConnectTimeout=20",
        remote_spec,
        str(local_parent),
    ]

    result = run_command(command)

    if result.returncode == 0:
        return True

    if "No such file or directory" in result.stderr:
        print(f"Remote path not found: {remote_folder}")
        return False

    print("scp failed.")
    print("Command:", " ".join(command))
    print("stderr:")
    print(result.stderr.strip())
    raise RuntimeError(f"scp failed with exit code {result.returncode}")


def shlex_quote(s: str) -> str:
    """Small wrapper so the rest of the code reads cleanly."""
    import shlex

    return shlex.quote(s)


# -----------------------------------------------------------------------------
# Job path helpers
# -----------------------------------------------------------------------------

def load_data_directory() -> Path:
    with DEFAULT_INPUT.open("r") as fp:
        input_json = json.load(fp)

    return Path(input_json["data_directory"]).expanduser()


def build_job_paths(
    data_directory: Path,
    job_set_name: str,
    attempt: int,
    m: int,
    n: int,
    temp: int,
) -> tuple[Path, Path, Path]:
    """
    Return local_job_folder, local_copy_parent, remote_job_folder.

    local_copy_parent is one level above T_<temp>, because recursive scp of a
    directory creates the final directory inside the destination.
    """
    jobfolder = local_job_parent(job_set_name)

    relative = (
        Path(jobfolder)
        / f"{job_set_name}_{attempt}"
        / f"M_{m}"
        / f"N_{n}"
        / f"T_{temp}"
    )

    remote_relative = (
        Path("jobs")
        / f"{job_set_name}_{attempt}"
        / f"M_{m}"
        / f"N_{n}"
        / f"T_{temp}"
    )

    local_job_folder = data_directory / relative
    local_copy_parent = local_job_folder.parent
    remote_job_folder = REMOTE_BASE_FOLDER / remote_relative

    return local_job_folder, local_copy_parent, remote_job_folder


def local_job_is_complete(local_job_folder: Path, n: int) -> bool:
    """Decide whether the job was already copied locally."""
    return (
        (local_job_folder / "timing.txt").exists()
        and (local_job_folder / "job_data.csv").exists()
        and u.find_max_index(str(local_job_folder) + os.sep) >= n
    )


def iter_jobs(config: TransferConfig, data_directory: Path) -> Iterable[tuple[Path, Path, Path, int]]:
    """Yield all requested jobs as local_folder, local_parent, remote_folder, n."""
    for job_set_name in config.job_set_names:
        for m in config.m_values:
            n = n_for_job(job_set_name, m, config.n_default, config.cbapa_c)
            for temp in config.temps:
                for attempt in config.attempts:
                    local_folder, local_parent, remote_folder = build_job_paths(
                        data_directory=data_directory,
                        job_set_name=job_set_name,
                        attempt=attempt,
                        m=m,
                        n=n,
                        temp=temp,
                    )
                    yield local_folder, local_parent, remote_folder, n


# -----------------------------------------------------------------------------
# Main transfer logic
# -----------------------------------------------------------------------------

def transfer_finished_jobs(config: TransferConfig) -> None:
    data_directory = load_data_directory()

    print(f"Local data directory: {data_directory}")
    print(f"Remote base folder:    {REMOTE_BASE_FOLDER}")
    print(f"Host alias:            {config.host_alias}")
    print(f"Dry run:               {config.dry_run}")
    print()

    copied = 0
    skipped_local = 0
    unfinished = 0

    # Fast sanity check: does the host alias work at all from Python?
    sanity = run_command(ssh_command(config, "echo SSH_OK"))
    if sanity.returncode != 0 or "SSH_OK" not in sanity.stdout:
        print("Could not connect with system ssh from inside Python.")
        print("stderr:")
        print(sanity.stderr.strip())
        raise RuntimeError("Initial ssh sanity check failed")

    for local_folder, local_parent, remote_folder, n in iter_jobs(config, data_directory):
        if str(local_folder) == str(remote_folder):
            raise RuntimeError("Remote and local folders are the same path; refusing to copy.")

        if local_job_is_complete(local_folder, n):
            print(f"Already copied: {local_folder}")
            skipped_local += 1
            continue

        remote_timing = remote_folder / "timing.txt"
        if not remote_file_exists(config, remote_timing):
            print(f"Not finished:    {remote_folder}")
            unfinished += 1
            continue

        print(f"Copying:        {remote_folder}")
        print(f"    to parent:  {local_parent}")

        if not config.dry_run:
            success = scp_get_folder(config, remote_folder, local_parent)
            if success:
                copied += 1
        else:
            copied += 1

    print()
    print("Summary")
    print(f"  copied/dry-run matches: {copied}")
    print(f"  already local:          {skipped_local}")
    print(f"  unfinished/missing:     {unfinished}")


def parse_args() -> TransferConfig:
    parser = argparse.ArgumentParser(description="Transfer finished Novus jobs to local SpaceLab data directory.")
    parser.add_argument("--host", default="Novus", help="Host alias from ~/.ssh/config. Default: Novus")
    parser.add_argument("--dry-run", action="store_true", help="Print what would be copied without copying.")
    parser.add_argument("--verbose-ssh", action="store_true", help="Add -v to ssh checks for debugging.")
    parser.add_argument(
        "--job-set",
        action="append",
        dest="job_sets",
        help="Job set name, e.g. BAPAWELD. Can be given multiple times.",
    )
    parser.add_argument("--temps", nargs="+", type=int, default=[1000], help="Temperatures to copy.")
    parser.add_argument(
        "--m-values",
        nargs="+",
        type=int,
        default=[3, 5, 10, 15, 20, 30, 50, 60, 100],
        help="M values to copy.",
    )
    parser.add_argument("--attempts", type=int, default=30, help="Number of attempts, starting from 0.")

    args = parser.parse_args()

    return TransferConfig(
        host_alias=args.host,
        job_set_names=tuple(args.job_sets) if args.job_sets else ("BAPAWELD",),
        attempts=range(args.attempts),
        m_values=tuple(args.m_values),
        temps=tuple(args.temps),
        dry_run=args.dry_run,
        verbose_ssh=args.verbose_ssh,
    )


def main() -> None:
    config = parse_args()
    transfer_finished_jobs(config)


if __name__ == "__main__":
    main()
