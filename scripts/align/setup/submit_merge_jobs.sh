#!/bin/bash

# Array of plate numbers and their corresponding S numbers
declare -A plate_s_numbers
plate_s_numbers[1]="1"
plate_s_numbers[2]="2"
plate_s_numbers[3]="4"
plate_s_numbers[4]="5"
plate_s_numbers[5]="3"

# Submit a job for each plate
for plate in "${!plate_s_numbers[@]}"; do
    sbatch scripts/cat_plates.sh "$plate" "${plate_s_numbers[$plate]}"
done