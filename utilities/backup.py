#!/usr/bin/env python3
import argparse
import os
import shutil
import subprocess
import sys
import time
from pathlib import Path
from datetime import datetime

def which_or_die(cmd):
    if shutil.which(cmd) is None:
        print(f"ERROR: '{cmd}' is not installed or not on PATH.", file=sys.stderr)
        sys.exit(1)

def run(cmd, check=True, capture=False, text=True):
    if capture:
        return subprocess.run(cmd, check=check, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=text)
    else:
        return subprocess.run(cmd, check=check)

def remote_exists(remote_name: str) -> bool:
    try:
        out = run(["rclone", "listremotes"], check=True, capture=True).stdout
        remotes = [r.strip().rstrip(":") for r in out.splitlines() if r.strip()]
        return remote_name in remotes
    except subprocess.CalledProcessError:
        return False

def ensure_dir(p: Path):
    p.mkdir(parents=True, exist_ok=True)

def mount_box(remote: str, mountpoint: Path, poll_interval="1m", dir_cache="12h", vfs_cache="writes", extra=[]):
    # Start rclone mount as a daemon
    cmd = [
        "rclone", "mount", f"{remote}:", str(mountpoint),
        "--daemon",
        "--vfs-cache-mode", vfs_cache,
        "--dir-cache-time", dir_cache,
        "--poll-interval", poll_interval,
        "--umask", "022",
        "--rc",  # allow control if needed
    ] + extra

    try:
        run(cmd, check=True)
    except subprocess.CalledProcessError as e:
        print(f"ERROR: rclone mount failed: {e}", file=sys.stderr)
        sys.exit(1)

    # Wait until the mount is ready
    for _ in range(60):  # up to ~30s
        # Prefer 'mountpoint -q' if available
        if shutil.which("mountpoint"):
            rc = subprocess.run(["mountpoint", "-q", str(mountpoint)]).returncode
            if rc == 0:
                return True
        # Fallback: try a simple list
        try:
            os.listdir(mountpoint)
            return True
        except Exception:
            pass
        time.sleep(0.5)

    print("ERROR: mount did not become ready in time.", file=sys.stderr)
    return False

def unmount(mountpoint: Path):
    # Try fusermount, then umount
    cmds = []
    if shutil.which("fusermount"):
        cmds.append(["fusermount", "-u", str(mountpoint)])
        cmds.append(["fusermount", "-uz", str(mountpoint)])  # lazy as fallback
    cmds.append(["umount", str(mountpoint)])
    for cmd in cmds:
        if subprocess.run(cmd).returncode == 0:
            return True
    return False

def build_rclone_copy_cmd(args, src: Path, remote_name: str, remote_subdir: str):
    base = [
        "rclone",
        # use 'copy' for no-deletes, 'sync' only in mirror mode
        "copy" if args.mode in ("new", "update") else "sync",
        str(src),
        f"{remote_name}:{remote_subdir}",
        "--stats=30s",
        "--progress" if args.progress else "--stats-one-line",
        "--transfers", str(args.transfers),
        "--checkers", str(args.checkers),
        "--retries", "5",
        "--low-level-retries", "10",
        "--retries-sleep", "10s",
        "--create-empty-src-dirs",
        "--metadata",  # keep metadata where supported
    ]

    # Mode flags
    if args.mode == "new":
        base.append("--ignore-existing")   # copy only if dest file doesn't exist
    elif args.mode == "update":
        # default 'copy' updates changed/new files; we keep it as-is,
        # but --update ensures we skip files newer on dest
        base.append("--update")

    # Mirror mode deletes — be explicit and cautious
    if args.mode == "mirror":
        base.extend(["--delete-during", "--fast-list"])

    # Safety / correctness
    if args.checksum:
        base.append("--checksum")  # slower; compare hashes instead of size+mtime
    if args.bwlimit:
        base.append(f"--bwlimit={args.bwlimit}")
    if args.dry_run:
        base.append("--dry-run")
    if args.fast_list:
        base.append("--fast-list")

    # Excludes
    for pat in (args.exclude or []):
        base.extend(["--exclude", pat])

    # Logging
    if args.log_file:
        base.extend(["--log-file", args.log_file, "--log-format", "date,time,microseconds", "--log-level", args.log_level])

    return base

def main():
    parser = argparse.ArgumentParser(
        description="Add-only backups to Box using rclone. Optionally mounts the Box remote for the duration."
    )
    parser.add_argument("--remote-name", default="box", help="Name of the rclone remote (default: box)")
    parser.add_argument("--mount", action="store_true", help="Mount the Box remote for this run")
    parser.add_argument("--mountpoint", default=str(Path.home() / "Box"), help="Where to mount the remote (default: ~/Box)")
    parser.add_argument("--keep-mounted", action="store_true", help="Keep the mount running after the backup")
    parser.add_argument("--local-subdir", required=True, help="Local source directory to back up")
    parser.add_argument("--remote-subdir", required=True, help="Remote subdirectory on Box (e.g. 'Backups/data')")
    parser.add_argument("--mode", choices=["new", "update", "mirror"], default="new",
                        help="new: copy files that do not exist on remote; "
                             "update: copy files when local is newer/different; "
                             "mirror: rclone sync (includes deletions!)")
    parser.add_argument("--exclude", action="append",
                        help="Exclude pattern (repeatable). Example: --exclude '.git/**' --exclude '*.tmp'")
    parser.add_argument("--bwlimit", type=str, default=None, help="Limit bandwidth (e.g., 20M, 800k)")
    parser.add_argument("--transfers", type=int, default=8, help="Concurrent file transfers")
    parser.add_argument("--checkers", type=int, default=8, help="Concurrent checking threads")
    parser.add_argument("--checksum", action="store_true", help="Use checksum comparison (slower, safest)")
    parser.add_argument("--dry-run", action="store_true", help="Show what would happen without changing anything")
    parser.add_argument("--progress", action="store_true", help="Show per-transfer progress")
    parser.add_argument("--fast-list", action="store_true", help="Use faster listing (more memory)")
    parser.add_argument("--log-file", default=None, help="rclone log file path")
    parser.add_argument("--log-level", default="INFO", choices=["DEBUG","INFO","NOTICE","ERROR"], help="rclone log level")
    args = parser.parse_args()

    # sanity
    which_or_die("rclone")

    # remote exists?
    if not remote_exists(args.remote_name):
        print(f"ERROR: rclone remote '{args.remote_name}' not found.\n"
              f"Create it with:  rclone config  (and add a remote named '{args.remote_name}' of type 'box')",
              file=sys.stderr)
        sys.exit(1)

    # paths
    src = Path(args.local_subdir).expanduser().resolve()
    if not src.is_dir():
        print(f"ERROR: Local path '{src}' does not exist or is not a directory.", file=sys.stderr)
        sys.exit(1)

    mount_started = False
    mountpoint = Path(args.mountpoint).expanduser().resolve()
    if args.mount:
        ensure_dir(mountpoint)
        print(f"[{datetime.now()}] Mounting {args.remote_name}: to {mountpoint} ...")
        if not mount_box(args.remote_name, mountpoint):
            print("ERROR: Failed to mount Box; aborting.", file=sys.stderr)
            sys.exit(1)
        mount_started = True
        print(f"[{datetime.now()}] Mount is ready.")

    # build and run rclone copy/sync
    cmd = build_rclone_copy_cmd(args, src, args.remote_name, args.remote_subdir.strip("/"))
    print(f"[{datetime.now()}] Running: {' '.join(cmd)}")
    rc = subprocess.run(cmd).returncode

    # unmount if we started it and user didn't request keep
    if mount_started and not args.keep-mounted:
        print(f"[{datetime.now()}] Unmounting {mountpoint} ...")
        ok = unmount(mountpoint)
        if not ok:
            print(f"WARNING: Unmount failed; you may need to run: fusermount -u {mountpoint}", file=sys.stderr)

    if rc != 0:
        print(f"rclone exited with code {rc}", file=sys.stderr)
        sys.exit(rc)
    print("[done]")

if __name__ == "__main__":
    main()
