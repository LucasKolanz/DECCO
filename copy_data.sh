#!/usr/bin/env bash

set -euo pipefail

usage() {
    cat << EOF
Usage:

  bash copy_jobs.sh \\
    --source 'SOURCE_GLOB' \\
    --destination 'DESTINATION_PATTERN' \\
    [--dry-run] [--ignore-existing] [--update]

Example:

  bash copy_jobs.sh \\
    --source '/mnt/drive/SpaceLab_data/jobsCosine/lognorm_*/N_300/T_1000' \\
    --destination '/mnt/drive/SpaceLab_data/jobs/BAPA_*/M_300/N_300/T_1000' \\
    --dry-run

Notes:

  The '*' part from the source is inserted into the '*' part of the destination.

  Example:
    source:      lognorm_6/N_300/T_1000
    destination: BAPA_6/M_300/N_300/T_1000
EOF
}

source_pattern=""
destination_pattern=""
dry_run=false
rsync_extra_opts=()

while [[ $# -gt 0 ]]; do
    case "$1" in
        --source)
            source_pattern="$2"
            shift 2
            ;;
        --destination|--dest)
            destination_pattern="$2"
            shift 2
            ;;
        --dry-run|--dryrun)
            dry_run=true
            shift
            ;;
        --ignore-existing)
            rsync_extra_opts+=(--ignore-existing)
            shift
            ;;
        --update)
            rsync_extra_opts+=(--update)
            shift
            ;;
        -h|--help)
            usage
            exit 0
            ;;
        *)
            echo "Unknown argument: $1"
            usage
            exit 1
            ;;
    esac
done

if [[ -z "$source_pattern" || -z "$destination_pattern" ]]; then
    echo "ERROR: --source and --destination are required."
    usage
    exit 1
fi

source_star_count="$(grep -o '\*' <<< "$source_pattern" | wc -l)"
destination_star_count="$(grep -o '\*' <<< "$destination_pattern" | wc -l)"

if [[ "$source_star_count" -ne 1 || "$destination_star_count" -ne 1 ]]; then
    echo "ERROR: this script currently expects exactly one '*' in --source and one '*' in --destination."
    echo "Source pattern:      $source_pattern"
    echo "Destination pattern: $destination_pattern"
    exit 1
fi

source_prefix="${source_pattern%%\**}"
source_suffix="${source_pattern#*\*}"

destination_prefix="${destination_pattern%%\**}"
destination_suffix="${destination_pattern#*\*}"

shopt -s nullglob

matches=( $source_pattern )

if [[ "${#matches[@]}" -eq 0 ]]; then
    echo "No matches found for:"
    echo "  $source_pattern"
    exit 0
fi

echo "Found ${#matches[@]} source match(es)."
echo

for src in "${matches[@]}"; do
    capture="$src"
    capture="${capture#"$source_prefix"}"
    capture="${capture%"$source_suffix"}"

    dst="${destination_prefix}${capture}${destination_suffix}"

    echo "============================================================"
    echo "SOURCE:      $src"
    echo "DESTINATION: $dst"
    echo "CAPTURE:     $capture"

    if [[ -d "$dst" ]]; then
        echo "DEST STATUS: exists"
    else
        if [[ "$dry_run" == true ]]; then
            echo "DEST STATUS: does not exist"
            echo "ACTION:      WOULD CREATE DESTINATION"
        else
            echo "DEST STATUS: does not exist"
            echo "ACTION:      CREATING DESTINATION"
            mkdir -p "$dst"
        fi
    fi

    echo "------------------------------------------------------------"

    if [[ "$dry_run" == true ]]; then
        rsync -ah  --dry-run \
            "${rsync_extra_opts[@]}" \
            "$src"/ "$dst"/
    else
        rsync -ah  \
            "${rsync_extra_opts[@]}" \
            "$src"/ "$dst"/
    fi

    echo
done

if [[ "$dry_run" == true ]]; then
    echo "Dry run only. No files or directories were changed."
fi