#!/bin/bash

# Define a cleanup function to run on interruption
#cleanup() {
#  echo "Script interrupted. Cleaning up..."
#  rm -r ./temp_canvas_images/
#  exit 1
#}

# Trap the SIGINT (Ctrl+C) signal and call the cleanup function
#trap cleanup SIGINT SIGTERM 

# Directory containing the ROOT files
root_dir="/afs/cern.ch/user/j/jsteenis/public/GEMS/gemDelays/results"

# Check if the directory exists
if [ ! -d "$root_dir" ]; then
    echo "Error: Directory $root_dir not found."
    exit 1
fi

# Iterate over each ROOT file in the directory
for root_file in "$root_dir"/finalFitInformation*.root; do
    if [ -f "$root_file" ]; then
        echo "Processing file: $root_file"
	python3 check_canvases.py ${root_file} ${root_file/.root/_check.pdf} 
	echo "Processing completed for $root_file"
    fi
done

for root_file in "$root_dir"/fitInformation*.root; do
    if [ -f "$root_file" ]; then
        echo "Processing file: $root_file"
	python3 check_canvases.py ${root_file} ${root_file/.root/_check.pdf} 
        echo "Processing completed for $root_file"
    fi
done

python3 check_means_canvases.py results/ results/check_means_final.pdf 
python3 check_means_canvases.py results/ results/check_means_initial.pdf "fitInformation*.root" 

echo "Processing all ROOT files in $root_dir completed."
