#This script runs the delay generation.
#To look under the hood of the classes, look at the delayClasses script
#Code by Jacob Steenis, 2024
from delayClasses_v2 import *
import glob

#files = ["/afs/cern.ch/user/j/jsteenis/public/GEMS/EfficiencyAnalyzer/results/delay_plots/gemPad_st1_Rneg1L2CH7_hist_chamberSeparated_fineYbinning.root"]
#files = glob.glob("./GEM_mcdonalds_data/*.root")
files = glob.glob("./GEM_mcdonalds_data/*M_14_L2.root")
#files = ["./results/delay_plots/GE11_P_10_L1_data.root"]

for i, input_file_name in enumerate(files):  
    print("Currently on file: ", input_file_name)  
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

    DG = delayGenerator(DR.histo, DR.histo_name, input_file_name, rebin_num=8, num_optimize_steps=5, reference_point=9)
    if DG.status == False:
        continue

    outfile = ROOT.TFile(f"GEM_delays/delays/{input_file_name.split('/')[-1].replace('.root','')}_delays.root", "RECREATE")
    outfile.cd()
    original_histo.Write()
    original_histo.ProfileX().Write(original_histo.GetName()+"_profileX")
    
    outfile.cd()
    float_applied_histo = DG.float_applied_histo
    int_applied_histo = DG.int_applied_histo
    int_differences = DG.int_differences
    gbt_applied_histo = DG.gbt_applied_histo
    gbt_differences = DG.gbt_differences
    
    float_applied_histo.Write()
    float_applied_histo.ProfileX().Write(float_applied_histo.GetName()+"_profileX")
    int_applied_histo.Write()
    int_applied_histo.ProfileX().Write(int_applied_histo.GetName()+"_profileX")
    int_differences.Write("integer_differences")
    gbt_applied_histo.Write()
    gbt_applied_histo.ProfileX().Write(gbt_applied_histo.GetName()+"_profileX")
    gbt_differences.Write("gbt_differences")
        
    outfile.Close()
    
    #DG.int_df.to_csv("results/delays.csv")
    if i==0:
        DG.group_df.to_csv("GEM_delays/delays/group_delays.csv", mode='w', columns=["padID", "fed", "amc", "oh", "vfat", "group", "bunchDelay"], index=False)
        DG.gbt_df.to_csv("GEM_delays/delays/gbt_delays.csv", mode='w', columns=["padID", "fed", "amc", "oh", "gbt", "gbtDelay"], index=False)

    else:
        DG.group_df.to_csv("GEM_delays/delays/group_delays.csv", mode='a', columns=["padID", "fed", "amc", "oh", "vfat", "group", "bunchDelay"], index=False, header=False)
        DG.gbt_df.to_csv("GEM_delays/delays/gbt_delays.csv", mode='a', columns=["padID", "fed", "amc", "oh", "gbt", "gbtDelay"], index=False, header=False)
