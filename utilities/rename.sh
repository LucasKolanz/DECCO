#!/bin/bash
# Rename Slurm jobs using a glob-with-captures pattern.
# Everything is literal EXCEPT '*' which:
#   • in PATTERN matches any sequence (including commas)
#   • in REPLACEMENT is replaced (left→right) by the corresponding capture
#
# Examples:
#   ./rename.sh 'relax,*' 'relax=1,*'
#   ./rename.sh '*,a=0,*' '*,a=1,*'
#
# NOTE: Always QUOTE the arguments so your shell doesn't expand the *.

set -u  # (no -e; we handle errors so we keep going)

if [ "$#" -lt 2 ] || [ "$#" -gt 3 ]; then
  echo "Usage: $0 'PATTERN' 'REPLACEMENT' [--dry-run]" >&2
  exit 1
fi

pattern="$1"
replacement="$2"
dry_run="${3-}"

# Count how many '*' (captures) are in the pattern
caps_only=${pattern//[^*]/}
num_caps=${#caps_only}

# Escape regex metachars in PATTERN, then turn '*' into '(.*)'
escape_regex() {
  # Escapes . [ ] { } ( ) + ? ^ $ | \ characters
  printf '%s' "$1" | sed -e 's/[.[\]{}()+?^$|\\]/\\&/g'
}
escaped_pat=$(escape_regex "$pattern")
pat_regex="^${escaped_pat//\*/(.*)}$"

updated=0; attempted=0

# Snapshot queue first so updates don't race the iterator
mapfile -t lines < <(squeue -u "${USER:-$(id -un)}" --noheader -o "%i|%j")

for line in "${lines[@]}"; do
  IFS='|' read -r jobid jobname <<<"$line"
  [[ -z "$jobid" || -z "$jobname" ]] && continue

  if [[ "$jobname" =~ $pat_regex ]]; then
    attempted=$((attempted+1))

    newname="$replacement"
    # Replace each '*' in REPLACEMENT with the corresponding capture
    for ((i=1; i<=num_caps; i++)); do
      cap="${BASH_REMATCH[$i]}"
      newname="${newname/\*/$cap}"
    done

    if [[ "$newname" != "$jobname" ]]; then
      echo "Renaming job $jobid: $jobname -> $newname"
      if [[ "$dry_run" == "--dry-run" ]]; then
        continue
      fi
      if scontrol update JobId="$jobid" JobName="$newname" 2>err.txt; then
        updated=$((updated+1))
      else
        echo "  WARN: scontrol failed for $jobid: $(<err.txt)" >&2
      fi
    fi
  fi
done

echo "Done. Updated $updated job(s) out of $attempted match(es)."
