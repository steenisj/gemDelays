import math
import pandas as pd
import numpy as np
import sys
import os
import re
import ROOT
import statistics

original_directory = os.getcwd()
sys.path.append('/afs/cern.ch/user/j/jsteenis/public/steenis-general-functions')
import general_functions as steen
sys.path.append(original_directory)

#For pulling the data to be processed
class dataRetriever():
    def __init__(self, file_path):
        self.file_path = file_path
        self.histo_name, self.histo = self.retriever()
        
    def retriever(self):
        input_file = ROOT.TFile.Open(self.file_path)
        for key in input_file.GetListOfKeys():
            if "histo_st" in str(key.GetName()):
                input_histo_name = str(key.GetName())
                input_histo = input_file.Get(str(key.GetName()))
                if int(input_histo.GetEntries()) == 0:
                    print("\033[91mThe file \033[0m", self.file_path.split("/")[-1], "\033[91mhas no entries!\033[0m")
                    return None, None
                else:
                    input_histo.SetDirectory(0)
                    return input_histo_name, input_histo
            else:
                raise ValueError("Was not able to locate the proper histo within .root file! :angerey_face:")
            
#This is initialized at the chamber level. Every chamber will have to re-initialize this class.
class delayGenerator():
    def __init__(self, histo, histo_name, filename, reference_point=10, rebin_num=8, division_width=64, num_optimize_steps=5, added_string="_intCorrectionApplied"):
        if histo==None or histo_name==None:
            self.status = False
        
        else:
            self.status = True
            self.histo = histo #Input histo
            self.histo_name = histo_name #Input histo name
            self.filename = filename
            self.reference_point = reference_point
            self.rebin_num = rebin_num
            self.division_width = division_width
            self.num_optimize_steps = num_optimize_steps
            self.added_string = added_string
            self.station, self.region, self.layer, self.chamber = self.gemPad_stringExtractor()
            self.all_df_float, self.expanded_differences_0, self.means_df, self.all_expanded_difference_hists = self.data_generator()            
            self.float_applied_histo = self.applier(self.all_df_float["0.0"]['delay'], self.histo, hist_string=self.added_string)
            
            self.int_df = self.int_optimizer()
            self.int_applied_histo = self.applier(self.int_df['bunchDelay'], self.histo, hist_string=self.added_string)
            self.int_differences = self.df_to_hist(self.int_df)
            
    def gemPad_stringExtractor(self):
        pattern = r"gemPad_st(\d+)_R([a-z]+)(\d+)L(\d+)CH(\d+)"
        match = re.match(pattern, self.filename.split("/")[-1])
            
        if match:
            station = int(match.group(1))
            region = match.group(2) + str(match.group(3))
            if "neg" in region:
                region = region.replace("neg","-")
                region = int(region)
            layer = int(match.group(4))
            chamber = int(match.group(5))
            return station, region, layer, chamber
        
        elif match is None:
            pattern = r"gemPad_st(\d+)_R(\d+)L(\d+)CH(\d+)"
            match = re.match(pattern, self.filename.split("/")[-1])
            
            if match is None:
                raise ValueError("Incorrect string format for the input!")
            
            station = int(match.group(1))
            region = str(match.group(2))
            layer = int(match.group(3))
            chamber = int(match.group(4))
            return station, region, layer, chamber   
        
    def calc_oh(self):
        oh = 2*((self.chamber-1)%6)
        return oh

    def calc_amc(self):
        amc = 2*int((self.chamber-1)/6.0) + 1
        return amc

    def calc_ifed(self):
        if self.station<0:
            return 1467
        elif self.station>0:
            return 1468
        else:
            raise ValueError("Somehow you got a station of 0 :(")

    def calc_vfat(self, df):
        vfat_values = []
        for padID in df['padID']:
            vfat = 8*int((padID%24)/192.0)+int(padID/64)
            vfat_values.append(vfat)
        return vfat_values
    
    def apply_info_to_df(self, df):
        df['oh'] = self.calc_oh()
        df['amc'] = self.calc_amc()
        df['vfat'] = self.calc_vfat(df)
    
    def grouper(self, hist):
        hist_copy = hist.Clone()
        return hist_copy.RebinX(self.rebin_num)
        
    def expander(self, hist):
        hist_copy = ROOT.TH1D(hist.GetName()+"_expanded", hist.GetName()+"_expanded", hist.GetNbinsX()*self.rebin_num, hist.GetXaxis().GetXmin(), hist.GetXaxis().GetXmax())
        
        #Expands the contents of one one bin in hist to be replicated for X num of bins in the copy        
        for i in range(1, (hist.GetNbinsX())*(self.rebin_num)+1):
            hist_copy.SetBinContent(i, hist.GetBinContent(math.ceil(i/(self.rebin_num))))
            #hist_copy.SetBinError(i, hist.GetBinError(math.ceil(i/(self.rebin_num))))     
        return hist_copy
    
    def data_generator(self):   
        all_expanded_difference_hists = {}
        all_delay_df = {}
        rebinned_hist = self.grouper(self.histo)
        fit_means_hist, fit_sigmas_hist = steen.fit_2d_histogram(rebinned_hist, "results/fitInformation_"+self.histo_name+".root", init_params=[0,8,1,0], param_limits={1:[0,15], 2:[0,4]})
        
        #Cycle through a range of reference points between self.reference_point and self.reference_point+1 to choose whichever removes the bias best        
        for i in np.linspace(0.0, 1, self.num_optimize_steps)[0:-1]:
            difference_hist = steen.compute_difference_histogram(fit_means_hist, self.reference_point+i, hist_name_str="_"+str(i))
            expanded_difference_hist = self.expander(difference_hist)
            all_expanded_difference_hists[i] = expanded_difference_hist
            all_delay_df[str(i)] = steen.histogram_to_df(expanded_difference_hist, "padID", "delay")
            self.apply_info_to_df(all_delay_df[str(i)])          
                    
        return all_delay_df, all_expanded_difference_hists[0], steen.histogram_to_df(fit_means_hist, "padID", "mean"), all_expanded_difference_hists
            
    def applier(self, correction_df, histo, hist_string=""):
        if len(correction_df)!=1536:
            print("Applied correction_df is not 1536 long!")
        applied_histo = ROOT.TH2D(histo.GetName()+hist_string, histo.GetName()+hist_string, histo.GetNbinsX(), histo.GetXaxis().GetXmin(), histo.GetXaxis().GetXmax(), histo.GetNbinsY(), histo.GetYaxis().GetXmin(), histo.GetYaxis().GetXmax())
        for n in range(1, histo.GetNbinsX()+1):
            for m in range(1, histo.GetNbinsY()+1):               
                #applied_histo.SetBinContent(n, m+int(self.expanded_differences_0.GetBinContent(n)*120), self.histo.GetBinContent(n, m))
                applied_histo.SetBinContent(n, m+int(correction_df.iloc[n-1]*120), histo.GetBinContent(n, m))
        return applied_histo
    
    def df_to_hist(self, df, histo_string=""):
        differences = ROOT.TH1D(self.histo.GetName()+"_intDifferences"+histo_string, self.histo.GetName()+"_intDifferences"+histo_string, self.histo.GetNbinsX(), self.histo.GetXaxis().GetXmin(), self.histo.GetXaxis().GetXmax())
        for i, pid in enumerate(df['padID']):
            differences.SetBinContent(differences.FindBin(pid), df['bunchDelay'][i])
        return differences
            
    def int_optimizer(self):
        all_delays_df_wInt = {}
        all_int_applied_histos = {}
        fit_stdev_values = {} #For storing the individual stdev values
        int_optimized_df = None
        
        for key in self.all_df_float.keys():
            DI = delayIntegerizer(self.all_df_float[key], self.histo, str(self.histo_name)+"_"+str(key))
            all_delays_df_wInt[key] = DI.delays
            all_int_applied_histos[key] = self.applier(DI.delays['bunchDelay'], self.histo, hist_string="_intApplied"+str(key).split(".")[-1])#DI.int_applied_histo
        
        #Cycle through all the vfats
        for vfat in np.arange(int(1536/self.division_width)):
            for i, diffs in enumerate(all_int_applied_histos.items()):
                temp_profileX = diffs[1].ProfileX()
                temp_binContent_values = [] #For storing the bin contents for this histo
                for j in range(int((self.division_width*vfat)+1), int((self.division_width*(vfat+1)))+1):
                    content = temp_profileX.GetBinContent(j)
                    temp_binContent_values.append(content)
                
                temp_stdev = statistics.stdev(temp_binContent_values)
                
                if vfat not in list(fit_stdev_values.keys()):
                    fit_stdev_values[vfat] = [temp_stdev]
                else:
                    fit_stdev_values[vfat].append(temp_stdev)
                    
        #Now that we have the correct reference_number_offset index, now we can collect the proper delays from the df
        print(fit_stdev_values)
        for key in fit_stdev_values.keys():  
            print(key)
            min_index = np.argmin(fit_stdev_values[key])
            min_offset_key = list(all_delays_df_wInt.keys())[min_index] 
            print(min_offset_key)
            
            if int_optimized_df is None:   
                int_optimized_df=all_delays_df_wInt[min_offset_key].iloc[int(key)*self.division_width:(int(key)+1)*self.division_width]
                
            else:
                int_optimized_df = pd.concat([int_optimized_df, all_delays_df_wInt[min_offset_key].iloc[int(key)*self.division_width:(int(key)+1)*self.division_width]], ignore_index=True)
        return int_optimized_df   
    
#For taking the float delays and converting them into the application-specific quantities, integerizing first.
class delayIntegerizer():
    def __init__(self, delay_df, histo, histo_name):                
        self.delays = delay_df
        self.histo = histo
        self.histo_name = histo_name
        self.delays_to_int() #Generates the integer corrections
        
    def delays_to_int(self):
        #self.delays['bunchDelay'] = self.delays['delay'].apply(lambda x: math.ceil(x))
        #self.delays['bunchDelay'] = self.delays['delay'].apply(lambda x: round(x, 0))
        temp_rounded_values = []
        for i, value in enumerate(self.delays['delay']):
            #if i<768:
            #    temp_rounded_values.append(round(value, 0))
            #else:
            temp_rounded_values.append(math.ceil(value))
                
        self.delays['bunchDelay'] = temp_rounded_values
        #print(self.delays['bunchDelay'], "\n", self.delays['delay'])
        self.delays['remainder'] = self.delays['bunchDelay'] - self.delays['delay']
