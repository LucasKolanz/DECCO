from pathlib import Path
from collections import defaultdict
import hashlib
import re
import csv
import json

# Project directory
relative_path = ""
project_path = Path(__file__).resolve().parent / relative_path

# Load default input
with open(project_path / "default_files/default_input.json", "r") as fp:
    input_json = json.load(fp)

# Base directory containing jobs/ and jobsCosine/
path = Path(input_json["data_directory"])

# Number of monomers in projectile
m = 20

source_pattern = "jobsCosine/lognorm_*/N_300/T_1000"
dest_pattern = f"jobs/DBAPA_*/M_{m}/N_20/T_1000"
file_pattern = f"{m}_*"


def file_hash(filename):
    """Calculate SHA-256 hash without loading the entire file into memory."""
    h = hashlib.sha256()

    with open(filename, "rb") as f:
        for chunk in iter(lambda: f.read(1024 * 1024), b""):
            h.update(chunk)

    return h.hexdigest()


def find_origins(path):
    # Only accept source folders named lognorm_<integer>
    source_folders = sorted(
        folder for folder in path.glob(source_pattern)
        if re.fullmatch(r"lognorm_\d+", folder.parent.parent.name)
    )

    dest_folders = sorted(path.glob(dest_pattern))

    print(f"Found {len(source_folders)} source folders")
    print(f"Found {len(dest_folders)} destination folders")

    # Index source files by (filename, hash)
    source_lookup = defaultdict(list)

    for folder in source_folders:
        print(f"Indexing {folder}")

        for file in folder.glob(file_pattern):
            if not file.is_file():
                continue

            key = (file.name, file_hash(file))
            source_lookup[key].append(folder)

    results = []

    for dest in dest_folders:
        matches = defaultdict(set)

        dest_files = sorted(
            file for file in dest.glob(file_pattern)
            if file.is_file()
        )

        for file in dest_files:
            key = (file.name, file_hash(file))

            for source in source_lookup.get(key, []):
                matches[source].add(file.name)

        # Find source folders with the most matching files
        best_count = max(
            (len(filenames) for filenames in matches.values()),
            default=0
        )

        best_sources = [
            source for source, filenames in matches.items()
            if len(filenames) == best_count
        ] if best_count > 0 else []

        if len(best_sources) == 1:
            status = "UNIQUE"
        elif len(best_sources) > 1:
            status = "AMBIGUOUS"
        else:
            status = "NOT FOUND"

        results.append({
            "destination": str(dest.relative_to(path)),
            "M": m,
            "sources": "; ".join(
                str(source.relative_to(path))
                for source in sorted(best_sources)
            ),
            "matched_files": best_count,
            "total_files": len(dest_files),
            "status": status,
        })

    return results


results = find_origins(path)

# Print results
for result in results:
    print(
        f"{result['destination']}\n"
        f"  Source: {result['sources'] or 'NOT FOUND'}\n"
        f"  Matching files: {result['matched_files']}/{result['total_files']}\n"
        f"  Status: {result['status']}\n"
    )

# Save results to CSV
output_file = path / f"projectile_origins_M_{m}.csv"

with open(output_file, "w", newline="") as f:
    writer = csv.DictWriter(
        f,
        fieldnames=[
            "destination",
            "M",
            "sources",
            "matched_files",
            "total_files",
            "status",
        ]
    )
    writer.writeheader()
    writer.writerows(results)

print(f"Results saved to {output_file}")
