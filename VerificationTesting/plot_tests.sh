#!/bin/bash

script_dir="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
path=$(jq -r '.data_directory' "$script_dir/../default_files/default_input.json")


python3 gen_weld_verification_plots.py  "$path"/testJobs/weldCollisionVerification_0 \
	--output weld_validation.png \
    --save-data weld_validation.csv

# python3 gen_big_box_verification_plots.py "$path"/testJobs/bigboxVerificationPeriodic_0/N_100/B_3.000e-03/T_3/ \
#     --output big_box_periodic_validation.png \
#     --save-data big_box_periodic_validation.csv



