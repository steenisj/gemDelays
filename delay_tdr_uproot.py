#This script generates the root histograms which will be used to make delays for the GEM chambers.
#Allows us to run the delay generation many times in a row with less computation time.

#The output is a histogram of expanded padID versus BX. 
#Note that the binning is in 1/120 of a BX to correspond to the increments at which a gbt can delay.
#Bins are centered at integers (or 1/120 of an integer)!

#Code by Jacob Steenis, 2024 

import uproot
import awkward as ak
import ROOT
from array import array
import numpy as np
import os
import glob

ROOT.gROOT.SetBatch(1)

def find_root_files(directory):
    root_files = []
    # Iterate over the directory and its subdirectories
    for root, dirs, files in os.walk(directory):
        # Filter files with ".root" extension
        root_files.extend(glob.glob(os.path.join(root, "*.root")))
    return root_files

def clusterIdExpander(firstId, numPads, etaPartition, **kwargs):
    output = [(firstId + x) + (8 - etaPartition)*192 for x in range(numPads)]
    #print("NumPads: ", numPads, "FirstID: ", firstId, "Output: ", output)
    return output

def clusterBXExpander(BX, numPads, **kwargs):
    output = [BX for x in range(numPads)]
    #print("NumPads: ", numPads, "PadBX: ", BX, "Output: ", output)
    return output

def shiftingBX(gemPadDigiCluster_PadBX, CSCConstants_LCT_CENTRAL_BX, tmbL1aWindowSize, gemPadDigiCluster_ClusterALCTMatchTime):
    gemBX = gemPadDigiCluster_PadBX + CSCConstants_LCT_CENTRAL_BX - int(tmbL1aWindowSize/2.0) - gemPadDigiCluster_ClusterALCTMatchTime 
    return gemBX

filelist = find_root_files("/eos/cms/store/group/dpg_gem/comm_gem/P5_Commissioning/cms-gem-automation/prod/prompt-v1/GEMCommonNTuples/385620")
tree_name = "muNtupleProducer/MuDPGTree"
filelist = [x+":"+tree_name for x in filelist]

variables = ["gemPadDigiCluster_layer", "gemPadDigiCluster_station", "gemPadDigiCluster_region", "gemPadDigiCluster_ClusterFirstPad", "gemPadDigiCluster_PadBX", "gemPadDigiCluster_PadClusterSize", "gemPadDigiCluster_etaPartition", "gemPadDigiCluster_chamber", "gemPadDigiCluster_ClusterALCTMatchTime"]

# Define the station, layer, and region values
stations = [1, 2]
layers = [1, 2]
regions = [-1, 1]
chambers = range(1,36)

CSCConstants_LCT_CENTRAL_BX, tmbL1aWindowSize = 8, 7

# Loop over station, layer, and region
for station in stations:
    for layer in layers:
        for region in regions:
            for chamber in chambers:
                #To prevent naming issues when accessing via root
                region_label = str(region)
                if region<0:
                    #region_label = str(region).replace("-","neg")
                    region_label = str(region).replace("-1","M")
                elif region>0:
                    region_label = str(region).replace("1","P")
                
                file_path = f"results/delay_plots/GE1{station}_{region_label}_{chamber}_L{layer}_data.root"
                hist_name = f"GE1{station}_{region_label}_{chamber}_L{layer}"

                #For checking to see if this has already been processed!
                if os.path.exists(file_path):
                    test_file = ROOT.TFile(file_path)
                    if not test_file.IsZombie():
                        if test_file.GetListOfKeys().Contains(hist_name):
                            if not test_file.Get(hist_name).IsZombie():
                                print(file_path, " exists and is not zombie, skipping it!")
                                continue

                # Create a TH2F histogram with ROOT
                #h2d = ROOT.TH2F(f"histo_st{station}_R{region_label}L{layer}CH{chamber}", "", int(1536), -0.5, 1535.5, 24*120, 0-0.004166666666666666, 24-0.004166666666666666)
                h2d = ROOT.TH2F(hist_name, "", int(1536), -0.5, 1535.5, 24*120, 0-0.004166666666666666, 24-0.004166666666666666)
                print(f"Station {station} Region {region_label} Layer {layer} Chamber {chamber}")
                
                # Loop over chunks of data using uproot.iterate
                counter = 0
                for data in uproot.iterate(filelist, variables, library='ak'):     
                    counter+=1           

                    location_CUT = (data['gemPadDigiCluster_region'] == region) & (data['gemPadDigiCluster_station'] == station) & (data['gemPadDigiCluster_layer'] == layer) & (data['gemPadDigiCluster_chamber'] == chamber)
                    data = data[location_CUT]
                    
                    gemBX = ak.Array([shiftingBX(gemPadDigiCluster_PadBX, CSCConstants_LCT_CENTRAL_BX, tmbL1aWindowSize, gemPadDigiCluster_ClusterALCTMatchTime) for gemPadDigiCluster_PadBX, gemPadDigiCluster_ClusterALCTMatchTime in zip(ak.flatten(data['gemPadDigiCluster_PadBX'],axis=None), ak.flatten(data['gemPadDigiCluster_ClusterALCTMatchTime'],axis=None))])
                    expandedClusters = ak.Array([clusterIdExpander(firstId, numPads, etaPartition) for firstId, numPads, etaPartition in zip(ak.flatten(data['gemPadDigiCluster_ClusterFirstPad'],axis=None), ak.flatten(data['gemPadDigiCluster_PadClusterSize'],axis=None), ak.flatten(data['gemPadDigiCluster_etaPartition']))])
                    expandedBX = ak.Array([clusterBXExpander(firstId, numPads) for firstId, numPads in zip(gemBX, ak.flatten(data['gemPadDigiCluster_PadClusterSize'],axis=None))])
                    
                    # Fill the histogram using FillN()
                    n_events = len(ak.flatten(expandedClusters, axis=None))
                    x_data = ak.flatten(expandedClusters, axis=None)
                    y_data = ak.flatten(expandedBX, axis=None)
                    weights = ak.ones_like(ak.flatten(expandedBX, axis=None))
                    
                    if n_events == 0:
                        print(f"No data for station{station}, layer{layer}, region{region_label}, chamber{chamber} for this file!")
                        continue
                    
                    h2d.FillN(n_events, array('d', x_data), array('d', y_data), array('d', weights))

                # Plot the histogram using ROOT
                c1 = ROOT.TCanvas("", "", 1200, 800)
                h2d.GetXaxis().SetTitle("Expanded Pad ID")
                h2d.GetYaxis().SetTitle("Bunch Crossing [Giovanni Eqn Applied]")
                h2d.Draw("COLZ")
                #h2d.SaveAs(f"results/delay_plots/GE1{station}_{region_label}_{chamber}_L{layer}_data.root")
                h2d.SaveAs(file_path)
