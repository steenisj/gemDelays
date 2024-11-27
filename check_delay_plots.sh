#!/bin/bash

# Initialize verbose mode to off
VERBOSE=0

# Function to print messages only in verbose mode
log_verbose() {
    if [ "$VERBOSE" -eq 1 ]; then
        echo "$@"
    fi
}

# Parse arguments
while [[ $# -gt 0 ]]; do
    case "$1" in
        -v|--verbose)
            VERBOSE=1
            ;;
        *)
            echo "Unknown option: $1"
            exit 1
            ;;
    esac
    shift
done

current_dir=$(pwd 2>/dev/null)
initial_root_dir="$current_dir"/GEM_delays/verification_plots/initial
final_root_dir="$current_dir"/GEM_delays/verification_plots/final
data_dir="$current_dir"/GEM_mcdonalds_data

# Check if the directory exists
if [ ! -d "$root_dir" ]; then
    echo "Error: Directory $root_dir not found."
    exit 1
fi

#For the verbose case, we output the individual canvas pdfs and the initial 2d distributions
if [[ "$VERBOSE" -eq 1 ]]; then
    echo
    log_verbose "Verbose mode is ON. Outputting the 2d initial data and the individual canvas fits"
    echo

    for root_file in "$final_root_dir"/finalFitInformation*.root; do
        if [ -f "$root_file" ]; then
            echo "Processing file: $root_file"
	    python3 check_canvases.py ${root_file} ${root_file/.root/_check.pdf} 
	    echo "Processing completed for $root_file"
        fi
    done

    echo

    for root_file in "$initial_root_dir"/fitInformation*.root; do
        if [ -f "$root_file" ]; then
            echo "Processing file: $root_file"
	    python3 check_canvases.py ${root_file} ${root_file/.root/_check.pdf} 
            echo "Processing completed for $root_file"
        fi
    done

    echo

    python3 check_initDistributions.py $data_dir "$initial_root_dir"/initial_mcdonalds_distributions.pdf

    echo
fi

python3 check_means_canvases.py $final_root_dir "$final_root_dir"/check_means_final.pdf 
echo
python3 check_means_canvases.py $initial_root_dir "$initial_root_dir"/check_means_initial.pdf "fitInformation*.root" 

echo "Processing all ROOT files in $root_dir completed."
