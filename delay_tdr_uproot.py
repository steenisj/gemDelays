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

variables = ["gemPadDigiCluster_layer", "gemPadDigiCluster_station", "gemPadDigiCluster_region", "gemPadDigiCluster_ClusterFirstPad", "gemPadDigiCluster_PadBX", "gemPadDigiCluster_PadClusterSize", "gemPadDigiCluster_etaPartition", "gemPadDigiCluster_chamber", "gemPadDigiCluster_ClusterALCTMatchTime", "mu_propagated_isME11", "mu_propagated_station", "mu_propagated_region", "mu_propagated_Outermost_z", "mu_propagated_etaP", "mu_propagated_pt", "mu_propagated_chamber", "mu_propagated_strip", "mu_propagated_layer"]

output_path = f"results/GEM_delay_data.root"

chamber_hists = {}

# Define the station, layer, and region values
stations = [1, 2]
layers = [1, 2]
regions = [-1, 1]
chambers = range(1,36+1)

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
                
                #file_path = f"results/delay_plots/GE{station}1_{region_label}_{chamber}_L{layer}_data.root"
                hist_name = f"GE{station}1_{region_label}_{chamber}_L{layer}"

                temp_hist = ROOT.TH2D(hist_name, hist_name, int(1536), -0.5, 1535.5, 24*120, 0-0.004166666666666666, 24-0.004166666666666666)
                temp_hist.GetXaxis().SetTitle("Expanded Pad ID")
                temp_hist.GetYaxis().SetTitle("Bunch Crossing [Giovanni Eqn Applied]")
                temp_hist.SetDirectory(0)

                chamber_hists[hist_name] = {"station": station, "region_label": region_label, "chamber": chamber, "layer": layer, "data": temp_hist} 

                #For checking to see if this has already been processed!
                '''if os.path.exists(file_path):
                    test_file = ROOT.TFile(file_path)
                    if not test_file.IsZombie():
                        if test_file.GetListOfKeys().Contains(hist_name):
                            if not test_file.Get(hist_name).IsZombie():
                                print(file_path, " exists and is not zombie, skipping it!")
                                continue'''

                # Create a TH2F histogram with ROOT
                #h2d = ROOT.TH2F(f"histo_st{station}_R{region_label}L{layer}CH{chamber}", "", int(1536), -0.5, 1535.5, 24*120, 0-0.004166666666666666, 24-0.004166666666666666)
                #h2d = ROOT.TH2F(hist_name, "", int(1536), -0.5, 1535.5, 24*120, 0-0.004166666666666666, 24-0.004166666666666666)
                #print(f"Station {station} Region {region_label} Layer {layer} Chamber {chamber}")
                



# Loop over chunks of data using uproot.iterate
counter = 0
counter_limit = None #10000
for data in uproot.iterate(filelist, variables, library='ak', step_size=1):     
    if counter_limit!=None:
        if counter>counter_limit:
            break
    counter+=1

    if (counter%1000==0):
        print("Working on event: ", counter)

    matchedChamber = []
    for muon in range(len(data['mu_propagated_isME11'])):
        
        if len(ak.flatten(data['mu_propagated_isME11'], axis=None)) == 0:
            continue
        
        if data['mu_propagated_isME11'][0][muon] and data['mu_propagated_station'][0][muon]==1 and data['mu_propagated_region'][0][muon]*data['mu_propagated_Outermost_z'][0][muon]>0 and data['mu_propagated_etaP'][0][muon]>=1 and data['mu_propagated_etaP'][0][muon]<=8 and data['mu_propagated_pt'][0][muon]>10:
           matchedChamber.append(data['mu_propagated_chamber'][0][muon]) 

    matchedChamber = list(dict.fromkeys(matchedChamber))

    if not matchedChamber:
        #print("skipping!")
        continue

    for muon in range(len(data['mu_propagated_isME11'])):
            muon_station = data['mu_propagated_station'][0][muon]
            muon_region = data['mu_propagated_region'][0][muon]
            muon_layer = data['mu_propagated_layer'][0][muon]
            muon_chamber = data['mu_propagated_chamber'][0][muon]

            if data['mu_propagated_isME11'][0][muon] and muon_station==1 and muon_region*data['mu_propagated_Outermost_z'][0][muon]>0 and data['mu_propagated_etaP'][0][muon]>=1 and data['mu_propagated_etaP'][0][muon]<=8 and data['mu_propagated_pt'][0][muon]>10:
                location_mask = (data['gemPadDigiCluster_region'] == muon_region) & ((data['gemPadDigiCluster_station'] == muon_station)) & (data['gemPadDigiCluster_layer'] == muon_layer) & (data['gemPadDigiCluster_chamber'] == muon_chamber) & (data['gemPadDigiCluster_etaPartition'] == data['mu_propagated_etaP'][0][muon])
                proximity_mask = (abs(data['gemPadDigiCluster_ClusterFirstPad'] - int(data['mu_propagated_strip'][0][muon]/2.0)) <= 5)
                
                for variable in variables:
                    if "gemPadDigiCluster" in variable:
                        data["temp_" + variable] = data[variable][location_mask & proximity_mask]

                gemBX = ak.Array([shiftingBX(gemPadDigiCluster_PadBX, CSCConstants_LCT_CENTRAL_BX, tmbL1aWindowSize, gemPadDigiCluster_ClusterALCTMatchTime) for gemPadDigiCluster_PadBX, gemPadDigiCluster_ClusterALCTMatchTime in zip(ak.flatten(data['temp_gemPadDigiCluster_PadBX'],axis=None), ak.flatten(data['temp_gemPadDigiCluster_ClusterALCTMatchTime'],axis=None))])
                expandedClusters = ak.Array([clusterIdExpander(firstId, numPads, etaPartition) for firstId, numPads, etaPartition in zip(ak.flatten(data['temp_gemPadDigiCluster_ClusterFirstPad'],axis=None), ak.flatten(data['temp_gemPadDigiCluster_PadClusterSize'],axis=None), ak.flatten(data['temp_gemPadDigiCluster_etaPartition']))])
                expandedBX = ak.Array([clusterBXExpander(firstId, numPads) for firstId, numPads in zip(gemBX, ak.flatten(data['temp_gemPadDigiCluster_PadClusterSize'],axis=None))])


                # Fill the histogram using FillN()
                n_events = len(ak.flatten(expandedClusters, axis=None))
                x_data = ak.flatten(expandedClusters, axis=None)
                y_data = ak.flatten(expandedBX, axis=None)
                weights = ak.ones_like(ak.flatten(expandedBX, axis=None))

                if n_events>0:
                    if muon_region == -1:
                        muon_region = "M"
                    elif muon_region == 1:
                        muon_region = "P"
                    name = f"GE{muon_station}1_{muon_region}_{muon_chamber}_L{muon_layer}"
                    chamber_hists[name]['data'].FillN(n_events, array('d', x_data), array('d', y_data), array('d', weights))

for chamber in chamber_hists.keys():
    chamber_hists[chamber]['data'].SaveAs("results/delay_plots/"+str(chamber)+".root")


