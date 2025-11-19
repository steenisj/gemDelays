#This script runs the delay generation.
#To look under the hood of the classes, look at the delayClasses script
#Code by Jacob Steenis, 2024/2025
from delayClasses_v2 import *
import glob

files = glob.glob("./GEM_mcdonalds_data/*.root")
#files = glob.glob("./GEM_mcdonalds_data/GE11_M_26_L2.root")
#files = glob.glob("./GEM_mcdonalds_data/GE11_P_34_L1.root")
#files = glob.glob("./GEM_mcdonalds_data/GE11_M_15_L1.root")

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

    DG = delayGenerator(DR.histo, DR.histo_name, input_file_name, rebin_num=8, num_optimize_steps=5, reference_point=7)
    
    if DG.status == False:
        continue

    outfile = ROOT.TFile(f"GEM_delays/delays/{input_file_name.split('/')[-1].replace('.root','')}_delays.root", "RECREATE")
    DG.histo.Write()
    DG.float_applied_histo.Write()
    DG.int_applied_histo.Write()
    DG.gbt_applied_histo.Write()

    DG.int_differences.Write("integer_differences")
    DG.gbt_differences.Write("gbt_differences")

    outfile.Close()
    
    if i==0:
        DG.group_df.to_csv(
                "GEM_delays/delays/group_delays.csv", 
                mode='w', 
                columns=["padID", "fed", "amc", "oh", "vfat", "group", "bunchDelay"], 
                index=False
            )

        DG.gbt_df.to_csv(
                "GEM_delays/delays/gbt_delays.csv", 
                mode='w', 
                columns=["padID", "fed", "amc", "oh", "gbt", "gbtDelay"], 
                index=False
            )

    else:
        DG.group_df.to_csv(
                "GEM_delays/delays/group_delays.csv", 
                mode='a', 
                columns=["padID", "fed", "amc", "oh", "vfat", "group", "bunchDelay"], 
                index=False, 
                header=False
            )

        DG.gbt_df.to_csv(
                "GEM_delays/delays/gbt_delays.csv", 
                mode='a', 
                columns=["padID", "fed", "amc", "oh", "gbt", "gbtDelay"], 
                index=False, 
                header=False
            )
    print("-------------------------------------------------------------------------------------------------------------\n\n\n")

print("\n-------------------------------------------------------------------------------------------------------------")
print("\033[1;35mPROCESS COMPLETED!\033[0m")
print("-------------------------------------------------------------------------------------------------------------\n")
