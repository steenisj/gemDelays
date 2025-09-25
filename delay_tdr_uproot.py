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
#filelist = find_root_files("/eos/cms/store/group/dpg_gem/comm_gem/P5_Commissioning/cms-gem-automation/prod/prompt-v1/GEMCommonNTuples/38562")
tree_name = "muNtupleProducer/MuDPGTree"
filelist = [x+":"+tree_name for x in filelist]

cluster_variables = ["gemPadDigiCluster_layer", "gemPadDigiCluster_station", "gemPadDigiCluster_region", "gemPadDigiCluster_ClusterFirstPad", "gemPadDigiCluster_PadBX", "gemPadDigiCluster_PadClusterSize", "gemPadDigiCluster_etaPartition", "gemPadDigiCluster_chamber", "gemPadDigiCluster_ClusterALCTMatchTime"]
muon_variables = ["mu_propagated_isME11", "mu_propagated_station", "mu_propagated_region", "mu_propagated_Outermost_z", "mu_propagated_etaP", "mu_propagated_pt", "mu_propagated_chamber", "mu_propagated_strip", "mu_propagated_layer"]
variables = cluster_variables + muon_variables

output_path = f"results/GEM_delay_data.root"
hot_output = ROOT.TFile("hotChannelsRemoved.root", "RECREATE")

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

                temp_hist = ROOT.TH2D(hist_name, hist_name, int(1536), -0.5, 1535.5, 24*120, 0-0.004166666666666666-0.5, 24-0.004166666666666666-0.5)
                temp_hist.GetXaxis().SetTitle("Expanded Pad ID")
                #temp_hist.GetYaxis().SetTitle("Bunch Crossing [Giovanni Eqn Applied]")
                temp_hist.GetYaxis().SetTitle("Time [bx]")
                temp_hist.SetDirectory(0)

                chamber_hists[hist_name] = {"station": station, "region_label": region_label, "chamber": chamber, "layer": layer, "data": temp_hist} 


# Loop over chunks of data using uproot.iterate
counter = 0
counter_limit = None #None if you want to run everything
for data in uproot.iterate(filelist, variables, library='ak', step_size=100000):     
    if counter_limit!=None:
        if counter>counter_limit:
            break
    counter+=1
    print("Working on event: ", counter*100000)


    muon_data = data[muon_variables]
    cluster_data = data[cluster_variables]
    #Check to make sure that we've got some sort of ME11 hit in the events
    ME11_mask = ak.any(muon_data['mu_propagated_isME11'], axis=-1)
    muon_data = muon_data[ME11_mask]
    cluster_data = cluster_data[ME11_mask]

    muon_location_mask = (muon_data['mu_propagated_isME11']) & (muon_data['mu_propagated_station']==1) & (muon_data['mu_propagated_region']*muon_data['mu_propagated_Outermost_z']>0) & (muon_data['mu_propagated_etaP']>=1) & (muon_data['mu_propagated_etaP']<=8) & (muon_data['mu_propagated_pt']>10)
    muon_data = muon_data[muon_location_mask]

    num_mask = ak.num(muon_data['mu_propagated_isME11'])>0
    muon_data = muon_data[num_mask]
    cluster_data = cluster_data[num_mask]

    chamber_combinations = ak.cartesian({"muon_chamber": muon_data['mu_propagated_chamber'], "cluster_chamber": cluster_data['gemPadDigiCluster_chamber']})
    layer_combinations = ak.cartesian({"muon_layer": muon_data['mu_propagated_layer'], "cluster_layer": cluster_data['gemPadDigiCluster_layer']})
    station_combinations = ak.cartesian({"muon_station": muon_data['mu_propagated_station'], "cluster_station": cluster_data['gemPadDigiCluster_station']})
    region_combinations = ak.cartesian({"muon_region": muon_data['mu_propagated_region'], "cluster_region": cluster_data['gemPadDigiCluster_region']})

    eta_combinations = ak.cartesian({"muon_eta": muon_data['mu_propagated_etaP'], "cluster_eta": cluster_data['gemPadDigiCluster_etaPartition']})
    proximity_combinations = ak.cartesian({"muon_prox": muon_data['mu_propagated_strip'], "cluster_prox": cluster_data['gemPadDigiCluster_ClusterFirstPad']})

    location_mask = (chamber_combinations['muon_chamber'] == chamber_combinations['cluster_chamber']) & (layer_combinations['muon_layer'] == layer_combinations['cluster_layer']) & (station_combinations['muon_station'] == station_combinations['cluster_station']) & (region_combinations['muon_region'] == region_combinations['cluster_region']) & (eta_combinations['muon_eta'] == eta_combinations['cluster_eta'])

    proximity_combinations['muon_prox'] = np.floor(proximity_combinations['muon_prox']/2.0) #To convert between strip and pad
    proximity_mask = (abs(proximity_combinations['cluster_prox'] - proximity_combinations['muon_prox']) <= 5)

    #Preparing the actual data
    BX_combinations = ak.cartesian({"muon_chamber": muon_data['mu_propagated_chamber'], "cluster_PadBX": cluster_data['gemPadDigiCluster_PadBX']})
    ALC_combinations = ak.cartesian({"muon_chamber": muon_data['mu_propagated_chamber'], "cluster_ALC": cluster_data['gemPadDigiCluster_ClusterALCTMatchTime']})
    size_combinations = ak.cartesian({"muon_chamber": muon_data['mu_propagated_chamber'], "cluster_size": cluster_data['gemPadDigiCluster_PadClusterSize']})
    firstPad_combinations = ak.cartesian({"muon_chamber": muon_data['mu_propagated_chamber'], "cluster_firstPad": cluster_data['gemPadDigiCluster_ClusterFirstPad']})

    BX_combinations = BX_combinations[location_mask & proximity_mask]
    ALC_combinations = ALC_combinations[location_mask & proximity_mask]
    size_combinations = size_combinations[location_mask & proximity_mask]
    eta_combinations = eta_combinations[location_mask & proximity_mask]
    firstPad_combinations = firstPad_combinations[location_mask & proximity_mask]

    chamber_combinations = chamber_combinations[location_mask & proximity_mask]
    layer_combinations = layer_combinations[location_mask & proximity_mask]
    station_combinations = station_combinations[location_mask & proximity_mask]
    region_combinations = region_combinations[location_mask & proximity_mask]

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
                    hist_station_mask = station_combinations['cluster_station'] == station
                    hist_layer_mask = layer_combinations['cluster_layer'] == layer
                    hist_region_mask = region_combinations['cluster_region'] == region
                    hist_chamber_mask = chamber_combinations['cluster_chamber'] == chamber
                    hist_location_mask = hist_station_mask & hist_layer_mask & hist_chamber_mask & hist_region_mask

                    gemBX = ak.Array([shiftingBX(gemPadDigiCluster_PadBX, CSCConstants_LCT_CENTRAL_BX, tmbL1aWindowSize, gemPadDigiCluster_ClusterALCTMatchTime) for gemPadDigiCluster_PadBX, gemPadDigiCluster_ClusterALCTMatchTime in zip(ak.flatten(BX_combinations['cluster_PadBX'][hist_location_mask],axis=None), ak.flatten(ALC_combinations['cluster_ALC'][hist_location_mask],axis=None))])
                    expandedClusters = ak.Array([clusterIdExpander(firstId, numPads, etaPartition) for firstId, numPads, etaPartition in zip(ak.flatten(firstPad_combinations['cluster_firstPad'][hist_location_mask],axis=None), ak.flatten(size_combinations['cluster_size'][hist_location_mask],axis=None), ak.flatten(eta_combinations['cluster_eta'][hist_location_mask]))])
                    expandedBX = ak.Array([clusterBXExpander(firstId, numPads) for firstId, numPads in zip(gemBX, ak.flatten(size_combinations['cluster_size'][hist_location_mask],axis=None))])

                    #hist_expandedClusters = expandedClusters[hist_location_mask]
                    #hist_expandedBX = expandedBX[hist_location_mask]

                    # Fill the histogram using FillN()
                    if len(expandedClusters) != len(expandedBX):
                        print("DIFFERENT LENGTHS OF CLUSTERS AND BXs")
                    n_events = len(ak.flatten(expandedClusters, axis=None))
                    x_data = ak.flatten(expandedClusters, axis=None)
                    y_data = ak.flatten(expandedBX, axis=None)
                    weights = ak.ones_like(ak.flatten(expandedBX, axis=None))

                    if n_events>0:
                        name = f"GE{station}1_{region_label}_{chamber}_L{layer}"
                        chamber_hists[name]['data'].FillN(n_events, array('d', x_data), array('d', y_data), array('d', weights))


print("REMOVING HOT CHANNELS")
for chamber in chamber_hists.keys():
    for i in range(chamber_hists[chamber]['data'].GetNbinsX()):
        threshold = 0
        for j in range(chamber_hists[chamber]['data'].GetNbinsY()):
            value = chamber_hists[chamber]['data'].GetBinContent(i,j)
            if value>threshold:
                threshold = value
        if threshold>10:
            threshold = 0.3*threshold
        else:
            threshold = 10e9


        high_bool = False
        low_bool = False
        override_bool = False
        if (chamber_hists[chamber]['data'].GetBinContent(i, chamber_hists[chamber]['data'].GetYaxis().FindBin(13)) >= threshold) or (chamber_hists[chamber]['data'].GetBinContent(i, chamber_hists[chamber]['data'].GetYaxis().FindBin(12)) >= threshold):
            high_bool = True
        if (chamber_hists[chamber]['data'].GetBinContent(i, chamber_hists[chamber]['data'].GetYaxis().FindBin(2)) >= threshold) or (chamber_hists[chamber]['data'].GetBinContent(i, chamber_hists[chamber]['data'].GetYaxis().FindBin(3)) >= threshold):
            low_bool = True
        if (chamber_hists[chamber]['data'].GetBinContent(i, chamber_hists[chamber]['data'].GetYaxis().FindBin(13)) > 2*threshold) and (chamber_hists[chamber]['data'].GetBinContent(i, chamber_hists[chamber]['data'].GetYaxis().FindBin(12)) > 2*threshold) and (chamber_hists[chamber]['data'].GetBinContent(i, chamber_hists[chamber]['data'].GetYaxis().FindBin(11)) > 2*threshold):
            override_bool = True

        if (high_bool and low_bool) or override_bool:
            print(f"THROWING AWAY chamber {chamber}, padID{i} SINCE IT'S HOT!")
            temp_hist = ROOT.TH1D(f"chamber{chamber}_padID{i}_removed", f"chamber{chamber}_padID{i}_removed", 24*120, 0-0.004166666666666666-0.5, 24-0.004166666666666666-0.5)
            for j in range(chamber_hists[chamber]['data'].GetNbinsY()):
                temp_hist.SetBinContent(j, chamber_hists[chamber]['data'].GetBinContent(i,j))
                #temp_hist.SetBinError(j, chamber_hists[chamber]['data'].GetBinError(i,j))
                chamber_hists[chamber]['data'].SetBinContent(i, j, 0)
                chamber_hists[chamber]['data'].SetBinError(i, j, 0)
            hot_output.cd()
            temp_hist.Write()

# Noisy pad filter
'''for chamber in chamber_hists.keys(): 
    max = 1
    for bx in range(chamber_hists[chamber]['data'].GetNbinsX()):
        for by in range(chamber_hists[chamber]['data'].GetNbinsY()):
            if chamber_hists[chamber]['data'].GetBinContent(bx, by) > max:
                max = chamber_hists[chamber]['data'].GetBinContent(bx, by)

        for i in [11,12,13,14,15,16]:
            if abs(max-chamber_hists[chamber]['data'].GetBinContent(bx, chamber_hists[chamber]['data'].GetYaxis().FindBin(i)))/max < 0.5:
                print(f"THROWING AWAY chamber {chamber}, padID{bx} SINCE IT'S HOT!")
                for by in range(chamber_hists[chamber]['data'].GetNbinsY()):
                    chamber_hists[chamber]['data'].SetBinContent(bx, by, 0)
                    chamber_hists[chamber]['data'].SetBinError(bx, by, 0)'''

hot_output.Close()
for chamber in chamber_hists.keys():
    chamber_hists[chamber]['data'].SaveAs("GEM_mcdonalds_data/"+str(chamber)+".root")


