#!/bin/bash

current_dir=$(pwd 2>/dev/null)
initial_root_dir="$current_dir"/GEM_delays/verification_plots/initial
final_root_dir="$current_dir"/GEM_delays/verification_plots/final

# Check if the directory exists
if [ ! -d "$root_dir" ]; then
    echo "Error: Directory $root_dir not found."
    exit 1
fi

# Iterate over each ROOT file in the directory
for root_file in "$final_root_dir"/finalFitInformation*.root; do
    if [ -f "$root_file" ]; then
        echo "Processing file: $root_file"
	python3 check_canvases.py ${root_file} ${root_file/.root/_check.pdf} 
	echo "Processing completed for $root_file"
    fi
done

for root_file in "$initial_root_dir"/fitInformation*.root; do
    if [ -f "$root_file" ]; then
        echo "Processing file: $root_file"
	python3 check_canvases.py ${root_file} ${root_file/.root/_check.pdf} 
        echo "Processing completed for $root_file"
    fi
done

python3 check_means_canvases.py $final_root_dir "$final_root_dir"/check_means_final.pdf 
python3 check_means_canvases.py $initial_root_dir "$initial_root_dir"/check_means_initial.pdf "fitInformation*.root" 

echo "Processing all ROOT files in $root_dir completed."
