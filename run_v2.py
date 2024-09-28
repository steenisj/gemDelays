from delayClasses_v2 import *

files = ["/afs/cern.ch/user/j/jsteenis/public/GEMS/EfficiencyAnalyzer/results/delay_plots/chamberSeparatedMcDonalds_rebinned/gemPad_st1_Rneg1L2CH6_hist_chamberSeparated_fineYbinning.root"]
num_iterations = 3

for input_file_name in files:    
    DR = dataRetriever(input_file_name)
    original_histo = DR.histo
    
    # rebin_num = number of pads you consider together 
    #       (8 pads = 1 group, whole group gets same correction)
    #
    # division_width means how many of these "rebinned pads" you minimize the post-correction st dev for 
    #       (here we do them all together, but you could say division_width=64 to consider the st dev at 
    #        a vfat level)
    #
    # num_optimize_step = number of increments between 0 and 1 you add to your reference number when 
    #       optimizing the differences from that reference number
    #
    # reference_point is the actual number you subtract the mean timing (per padID) to get estimates for
    #       the delays needed
    
    outfile = ROOT.TFile(f"results/{input_file_name.split('/')[-1].replace('.root','')}_delays.root", "RECREATE")
    outfile.cd()
    original_histo.Write()
    
    for i in range(0,num_iterations):        
        if i==0:    
            DG = delayGenerator(DR.histo, DR.histo_name, input_file_name, rebin_num=8, division_width=1536, num_optimize_steps=10, reference_point=10)    
            if not DG.status:
                continue
            
        else:
            prev_float_applied_histo = float_applied_histo
            prev_int_applied_histo = int_applied_histo
            prev_int_differences = int_differences
            
            DG = delayGenerator(prev_int_applied_histo, DR.histo_name, input_file_name, rebin_num=8, division_width=1536, num_optimize_steps=10, reference_point=10+i, added_string="")    
            if not DG.status:
                continue
            
        outfile.cd()
        float_applied_histo = DG.float_applied_histo
        int_applied_histo = DG.int_applied_histo
        int_differences = DG.int_differences
    
        float_applied_histo.Write(int_applied_histo.GetName()+str(i))
        int_applied_histo.Write(int_applied_histo.GetName()+str(i))
        int_differences.Write("integer_differences"+str(i))
        
    outfile.Close()
