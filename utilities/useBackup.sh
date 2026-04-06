#!/bin/bash

python3 backup.py \
--remote-name box \
--remote-subdir SpaceLab_data/ \
--local-subdir /mnt/49f170a6-c9bd-4bab-8e52-05b43b248577/SpaceLab_data/ \
--log-file /mnt/49f170a6-c9bd-4bab-8e52-05b43b248577/SpaceLab/utilities/logs/rcloneBoxlogs.log \
--mode update \

# python3 backup.py \
# --remote-name onedrive \
# --remote-subdir SpaceLab_data/ \
# --local-subdir /mnt/49f170a6-c9bd-4bab-8e52-05b43b248577/SpaceLab_data/ \
# --log-file /mnt/49f170a6-c9bd-4bab-8e52-05b43b248577/SpaceLab/utilities/logs/rcloneOneDrivelogs.log \
# --mode update \