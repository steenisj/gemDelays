This repo is designed to generate the GEM delays based on un-delayed trigger data. The general strategy for generating the delays comes in two parts: the group-level delays and the gbt-level delays. 

For the group-level (integer BX) delays, the approach is to choose a reference number---typically 7---from which we subtract the fitted mean group timing which we round to an integer. Because the intgerization of this difference is relatively arbitrary and can--in some cases--lead to worse performance than with no delay, we shift the reference number by increments of 0.2 BXs minimize the standard deviation between fitted group means on a per-chamber basis. In other words, for a reference number of 7, we check the resulting mcdonalds plots for 7.0, 7.2, 7.4, 7.6, and 7.8 and choose the reference number which leads to the smallest standard deviation of the group timings. We then integerize the differences from THIS reference number.

Now, this leads to the center of the distribution being offset by some n*0.2 amount. To fix this, we use the gbt delays (which pull the timing in the opposite direction) to offset this amount, centering the mean-timing distributions back to the original reference number --- in this example 7.

This process is run in two parts:
1) Execute: <br><br>```python3 generate_mcdonalds_plots.py -r [run number] -l [optional data limit: (n+1)*100,000 events within run]```<br><br> to generate the "McDonalds" 2d plots of padID vs bx. The outputs will be generated in a folder called GEM_mcdonalds_data.<br>
2) Execute: <br><br>```python3 run.py -r [run number]```<br><br> to generate text (and root) outputs for the delays. The outputs will be in a newly generated folder called GEM_delays.

Summary of outputs:
1) For the first part of the processing, files appear in the run[run number] folder within GEM_mcdonalds_data directory. The files look like, for example, GE11_P_9_L2.root. The only object within is a TH2D named, for example, GE11_P_9_L2. This is the McDonalds plot for chamber GE11_P_9_L2.
   
2) There are a lot more outputs for the second stage of the processing (as validatation checks and plots). Generally, they are piled into various subfolders within the GEM_delays directory.<br><br>
  a) delays/gbt_delays.csv --- These are the gbt level delays as found in the run.py script.<br>
  b) delays/group_delays.csv --- These are the delays for groups of 8 pads on the chambers.<br>
  c) delays/GE11_M_29_L2_delays.root (for example) --- These plots contain the vast majority of the correction information. Inside, we have the following objects. For TH2D objects, we have **GE11_M_29_L2** (the original data histogram, binned in groups); **GE11_M_29_L2_floatCorrectionApplied** (which is an ideal case where we could perfectly delay the data (relative to the reference value) by floats in the groups); **GE11_M_29_L2_intCorrectionApplied** (which is the data after the group delays are applied); **GE11_M_29_L2_intCorrectionApplied_gbtCorrectionApplied** (which is the data after the gbt delays are applied). For TH1D objects, we have **integer_differences** (the group delays); and **gbt_differences** (the gbt delays). <br>
  d) verification_plots/inital or verification_plots/final directories --- These directories will be filled with all the fit information from each chamber so that you can verify cases where you suspect the initial (pre-correction) or final (post-correction) fits are wonky.<br>
  e) verification_plots/intermediate directory --- Generally, these plots contain information about the fits for the optimization step of the delay correction. The index at the end of the name corresponds to the index within [0,0.2,0.4,0.6,0.8] as the additions onto the reference number.<br>
