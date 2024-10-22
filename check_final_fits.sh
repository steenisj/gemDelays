# Directory containing the ROOT files
root_dir="/afs/cern.ch/user/j/jsteenis/public/GEMS/gemDelays/results"

# Check if the directory exists
if [ ! -d "$root_dir" ]; then
    echo "Error: Directory $root_dir not found."
    exit 1
fi

# Iterate over each ROOT file in the directory
for root_file in "$root_dir"/*finalFitInformation_histo*; do
#for root_file in "$root_dir"/*_peakFit2d.root; do
    # Check if there are any matching files
    if [ -f "$root_file" ]; then
        echo "Processing file: $root_file"
        python3 check_canvases.py ${root_file} ${root_file/.root/_check.pdf}
        echo "Processing completed for $root_file"
    fi
done

echo "Processing all ROOT files in $root_dir completed."
