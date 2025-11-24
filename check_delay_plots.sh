#!/bin/bash

# Initialize verbose mode to off
VERBOSE=0
CHECK_INITIAL_ONLY=0
RUN_NUMBER=-1

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
	--check_initial)
	    CHECK_INITIAL_ONLY=1
	    ;;
	-r|--run)
            if [[ -n "$2" && ! "$2" =~ ^- ]]; then
                RUN_NUMBER="$2"
                shift
            else
                echo "Error: --run requires a numeric argument."
                return 1
            fi
            ;;

        *)
            echo "Unknown option: $1"
            return 1
            ;;
    esac
    shift
done

current_dir=$(pwd 2>/dev/null)

# Determine subfolder name
if [[ "$RUN_NUMBER" != "-1" ]]; then
    RUN_SUBFOLDER="run${RUN_NUMBER}"
else
    RUN_SUBFOLDER="default"
fi

# Construct paths
initial_root_dir="$current_dir/GEM_delays/${RUN_SUBFOLDER}/verification_plots/initial"
final_root_dir="$current_dir/GEM_delays/${RUN_SUBFOLDER}/verification_plots/final"
delays_dir="$current_dir/GEM_delays/${RUN_SUBFOLDER}/delays"
mcdonalds_dir="$current_dir/GEM_mcdonalds_data/${RUN_SUBFOLDER}"

if [[ "$CHECK_INITIAL_ONLY" -eq 1 ]]; then
    echo
    echo "Running only initial 2D distribution check..."
    python3 check_2d_distributions.py "$mcdonalds_dir" "$mcdonalds_dir"/initial_mcdonalds_distributions.pdf "GE*.root" ""
    echo
    return 0
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

    python3 check_2d_distributions.py $delays_dir "$initial_root_dir"/initial_mcdonalds_distributions.pdf
    echo
    python3 check_2d_distributions.py $delays_dir "$final_root_dir"/final_mcdonalds_distributions.pdf "GE*_delays.root" "_intCorrectionApplied_gbtCorrectionApplied"
    echo
fi

python3 check_means_canvases.py $final_root_dir "$final_root_dir"/check_means_final.pdf "finalFitInformation*.root" 
echo
python3 check_means_canvases.py $initial_root_dir "$initial_root_dir"/check_means_initial.pdf "fitInformation*.root"
echo
python3 check_means_canvases.py $initial_root_dir "$initial_root_dir"/check_sigmas_initial.pdf "fitInformation*.root" "fit_sigmas_hist" [0,10] 
echo
python3 check_means_canvases.py $initial_root_dir "$initial_root_dir"/check_amplitudes_initial.pdf "fitInformation*.root" "fit_amplitudes_hist" [0,500]
echo
python3 check_means_canvases.py $initial_root_dir "$initial_root_dir"/check_backgrounds_initial.pdf "fitInformation*.root" "fit_backgrounds_hist" [0,500]
echo

#echo "Processing all ROOT files in $root_dir completed."
