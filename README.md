This repo is designed to generate the GEM delays based on un-delayed trigger data.

It is run in two parts:
1) Execute ```python3 generate_mcdonalds_plots.py -r [run number] -l [optional data limit: (n+1)*100,000 events within run]``` to generate the "McDonalds" 2d plots of padID vs bx. The outputs will be generated in a folder called GEM_mcdonalds_data.
2) Execute ```python3 run.py -r [run number]``` to generate text (and root) outputs for the delays. The outputs will be in a newly generated folder called GEM_delays.
