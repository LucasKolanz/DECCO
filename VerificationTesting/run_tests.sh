#!/bin/bash

script_dir="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"



#welding verification tests
python3 "$script_dir"/make_weld_verification.py -r


#big box verification tests
#big box periodic test
# python3 "$script_dir"/make_big_box_periodic_verification.py -r
