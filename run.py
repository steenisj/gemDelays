#This script runs the delay generation.
#To look under the hood of the classes, look at the delayClasses script
#Code by Jacob Steenis, 2024/2025
from delayClasses import *
import glob
import argparse
from setup import ensure_folders_exist

def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('-r', '--run', type=int, default=None, help='Process a single run number')
    args = parser.parse_args()
    return args

if __name__ == '__main__':
    args = parse_args()
    ensure_folders_exist("temp_canvas_images", [])

    subfolders = [
                    "delays",
                    "verification_plots",
                    "verification_plots/final",
                    "verification_plots/initial",
                    "verification_plots/intermediate"
                ]

    if args.run is not None:
        ensure_folders_exist(f"GEM_delays/run{args.run}", subfolders)
        files = glob.glob(f"./GEM_mcdonalds_data/run{args.run}/*.root") #pull all files in the GEM_mcdonalds file for a specific run

    else:
        ensure_folders_exist("GEM_delays/default", subfolders)
        files = glob.glob("./GEM_mcdonalds_data/default/*.root") #pull all files in the GEM_mcdonalds file

    all_hot_channels = {}
    for i, input_file_name in enumerate(files): 
        print("-------------------------------------------------------------------------------------------------------------")
        print("\n\033[1;32mCurrently on file: \033[0m", input_file_name)
        DR = dataRetriever(input_file_name)
        original_histo = DR.histo
    
        # rebin_num = number of pads you consider together 
        #       (8 pads = 1 group, whole group gets same correction)
        #
        # num_optimize_step = number of increments between 0 and 1 you add to your reference number when 
        #       optimizing the differences from that reference number
        #
        # reference_point is the actual number you subtract the mean timing (per padID) to get estimates for
        #       the delays needed

        DG = delayGenerator(DR.histo, DR.histo_name, input_file_name, rebin_num=8, num_optimize_steps=5, reference_point=7, SPECIFY_RUN=args.run)
    
        if DG.status == False:
            continue

        outfile = ROOT.TFile(f"{DG.baseName}/delays/{input_file_name.split('/')[-1].replace('.root','')}_delays.root", "RECREATE")
        DG.histo.Write()

        for chamber, channels in DG.hotChannels.items():
            if chamber not in all_hot_channels:
                all_hot_channels[chamber] = set()

            all_hot_channels[chamber].update(channels)

        DG.float_applied_histo.Write()
        DG.int_applied_histo.Write()
        DG.gbt_applied_histo.Write()

        DG.int_differences.Write("integer_differences")
        DG.gbt_differences.Write("gbt_differences")

        outfile.Close()
    
        if i==0:
            DG.group_df.to_csv(
                    f"{DG.baseName}/delays/group_delays.csv", 
                    mode='w', 
                    columns=["padID", "fed", "amc", "oh", "vfat", "group", "bunchDelay"], 
                    index=False
                )

            DG.gbt_df.to_csv(
                    f"{DG.baseName}/delays/gbt_delays.csv", 
                    mode='w', 
                    columns=["padID", "fed", "amc", "oh", "gbt", "gbtDelay"], 
                    index=False
                )

        else:
            DG.group_df.to_csv(
                    f"{DG.baseName}/delays/group_delays.csv", 
                    mode='a', 
                    columns=["padID", "fed", "amc", "oh", "vfat", "group", "bunchDelay"], 
                    index=False, 
                    header=False
                )

            DG.gbt_df.to_csv(
                    f"{DG.baseName}/delays/gbt_delays.csv", 
                    mode='a', 
                    columns=["padID", "fed", "amc", "oh", "gbt", "gbtDelay"], 
                    index=False, 
                    header=False
                )
        print("-------------------------------------------------------------------------------------------------------------\n")

        with open(f"{DG.baseName}/all_hot_channels.txt", "w") as f:
            for chamber, channels in all_hot_channels.items():
                f.write(f"{chamber}: {sorted(channels)}\n")

    print("\n-------------------------------------------------------------------------------------------------------------")
    print("\033[1;35mPROCESS COMPLETED!\033[0m")
    print("-------------------------------------------------------------------------------------------------------------\n")
