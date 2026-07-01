#!/usr/bin/env bash
set -euo pipefail

usage() {
    cat <<EOF
Usage:
  $(basename "$0") PATTERN [--yes|-y] [--ignore-case|-i]

Cancels your Slurm jobs whose job name contains PATTERN.

By default, this does a dry run and only prints matching jobs.
Use --yes or -y to actually cancel them.

Examples:
  $(basename "$0") DBAPA
  $(basename "$0") relax --yes
  $(basename "$0") relax --ignore-case --yes
EOF
}

if [[ $# -lt 1 ]]; then
    usage
    exit 1
fi

pattern=""
do_cancel=false
ignore_case=false

while [[ $# -gt 0 ]]; do
    case "$1" in
        --yes|-y)
            do_cancel=true
            shift
            ;;
        --ignore-case|-i)
            ignore_case=true
            shift
            ;;
        --help|-h)
            usage
            exit 0
            ;;
        -*)
            echo "Unknown option: $1" >&2
            usage
            exit 1
            ;;
        *)
            if [[ -z "$pattern" ]]; then
                pattern="$1"
                shift
            else
                echo "Unexpected extra argument: $1" >&2
                usage
                exit 1
            fi
            ;;
    esac
done

if [[ -z "$pattern" ]]; then
    echo "Error: missing PATTERN" >&2
    usage
    exit 1
fi

if "$ignore_case"; then
    matches=$(
        squeue -u "$USER" --noheader -o "%i|%j" \
            | awk -F'|' -v pat="$pattern" '
                index(tolower($2), tolower(pat)) {print $1 "|" $2}
              '
    )
else
    matches=$(
        squeue -u "$USER" --noheader -o "%i|%j" \
            | awk -F'|' -v pat="$pattern" '
                index($2, pat) {print $1 "|" $2}
              '
    )
fi

if [[ -z "$matches" ]]; then
    echo "No jobs found with pattern: $pattern"
    exit 0
fi

echo "Matching jobs:"
echo "$matches" | awk -F'|' '{printf "  %s  %s\n", $1, $2}'

if ! "$do_cancel"; then
    echo
    echo "Dry run only. Re-run with --yes to cancel these jobs."
    exit 0
fi

echo
echo "Cancelling matching jobs..."
echo "$matches" | awk -F'|' '{print $1}' | xargs -r scancel
echo "Done."