#These are the classes that are called when generating the GEM delays.
#Note that all the work is done within the __init__ function of the delayGenerator ... look there first!
#Code by Jacob Steenis, 2024
import math
import pandas as pd
import numpy as np
import sys
import os
import re
import ROOT
import statistics

#This is initialized at the chamber level. Every chamber will have to re-initialize this class.
class delayGenerator():
    def __init__(self, histo, histo_name, filename, reference_point=10, rebin_num=8, num_optimize_steps=5):
        if histo==None or histo_name==None:
            self.status = False
        
        else:
            self.status = True
            
            #Getting parameters
            self.general = generalFunctions()
            self.histo = histo #Input histo
            self.histo_name = histo_name #Input histo name
            self.filename = filename
            self.reference_point = reference_point
            self.rebin_num = rebin_num
            self.num_optimize_steps = num_optimize_steps
            
            #Extracting the data and generating what the "ideal delays" would be
            self.station, self.region, self.layer, self.chamber = self.gemPad_stringExtractor()
            self.all_df_float, self.expanded_differences_0, self.means_df, self.all_expanded_difference_hists = self.data_generator()            
            self.float_applied_histo = self.applier(self.all_df_float["0.0"]["idealDelay"], self.histo, hist_string="_floatCorrectionApplied")
            
            #Getting the integer delays for the groups
            self.int_df, self.min_reference_point = self.int_optimizer()
            self.int_applied_histo = self.applier(self.int_df["bunchDelay"], self.histo, hist_string="_intCorrectionApplied")
            self.int_differences = self.df_to_hist(self.int_df["bunchDelay"], self.int_df["padID"], histo_string="_intDifferences")
            
            #Adding in the gbt delays
            self.final_df = self.calc_gbt_delay() #Adds gbt delays to the wf
            self.gbt_applied_histo = self.gbt_applier(self.final_df['gbtDelay'], self.int_applied_histo, hist_string="_gbtCorrectionApplied")
            self.gbt_differences = self.df_to_hist(self.final_df["gbtDelay"], self.final_df["padID"], histo_string="_gbtDifferences")
            
            #Formats the output histos
            self.format_histos()

            #This is just for checking the effect of the delays
            self.final_means, self.final_sigmas = self.general.fit_2d_histogram(self.gbt_applied_histo, 
                                                                                output_file="results/finalFitInformation_"+self.gbt_applied_histo.GetName()+".root", 
                                                                                init_params=[0,9,2,0], 
                                                                                param_limits={1:[4,15], 2:[0,4]}, 
                                                                                fit_range=[5,13], 
                                                                                pol0_from_back=14
                                                                                )
            
            #final = with all information
            self.final_df["padID"] = self.final_df["padID"].astype(int)
            self.final_df["bunchDelay"] = self.final_df["bunchDelay"].astype(int)
            
            #For formatting to be input into the electronics (reducing the final_df)
            self.group_df = self.df_reducer_group(self.final_df)
            self.gbt_df = self.df_reducer_gbt(self.final_df)
    
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
        if self.layer== 2:
            oh += 1
        return oh

    def calc_amc(self):
        amc = 2*int((self.chamber-1)/6.0) + 1
        return amc

    def calc_ifed(self):
        if self.region<0:
            return 1467
        elif self.region>0:
            return 1468
        else:
            raise ValueError("Somehow you got a station of 0 :(")

    def calc_vfat(self, df):
        vfat_values = []
        for padID in df['padID']:
            eta = int(padID/192)
            phi = int((padID%192)/64)
            #vfat = 8*int((padID%24)/192.0)+int(padID/64)
            vfat = eta+phi*8 
            vfat_values.append(vfat)
        return vfat_values

    def calc_gbt(self, df):
        gbt0 = [7,12,13,14,15,23]
        gbt1 = [0,1,2,3,4,5,6,8,16]
        gbt2 = [9,10,11,17,18,19,20,21,22]
        
        gbt_values = []
        
        for vfat in df['vfat']:
            if int(vfat) in gbt0:
                temp_gbt = 0
            elif int(vfat) in gbt1:
                temp_gbt = 1
            elif int(vfat) in gbt2:
                temp_gbt = 2
            else:
                temp_gbt = -9999999999
            gbt_values.append(temp_gbt)
            
        return gbt_values
            
    def calc_group(self, df):
        group_values = []
        for padID in df['padID']:
            group = int(padID/self.rebin_num)%8
            group_values.append(group)
        return group_values
    
    def apply_info_to_df(self, df):
        df['oh'] = self.calc_oh()
        df['amc'] = self.calc_amc()
        df['fed'] = self.calc_ifed()
        df['vfat'] = self.calc_vfat(df)
        df['gbt'] = self.calc_gbt(df)
        df['group'] = self.calc_group(df)
    
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
        fit_means_hist, fit_sigmas_hist = self.general.fit_2d_histogram(rebinned_hist, output_file="results/fitInformation_"+self.histo_name+".root", init_params=[0,7,2,0], param_limits={1:[2,13], 2:[0,4]}, fit_range=[3,11], pol0_from_back=12)
        
        #Cycle through a range of reference points between self.reference_point and self.reference_point+1 to choose whichever removes the bias best        
        for i in np.linspace(0.0, 1.0, self.num_optimize_steps+1)[0:-1]:
            difference_hist = self.general.compute_difference_histogram(fit_means_hist, self.reference_point+i, hist_name_str="_"+str(i))
            expanded_difference_hist = self.expander(difference_hist)
            all_expanded_difference_hists[i] = expanded_difference_hist
            all_delay_df[str(i)] = self.general.histogram_to_df(expanded_difference_hist, "padID", "idealDelay")
            self.apply_info_to_df(all_delay_df[str(i)])          
                    
        return all_delay_df, all_expanded_difference_hists[0], self.general.histogram_to_df(fit_means_hist, "padID", "mean"), all_expanded_difference_hists
            
    def applier(self, correction_df, histo, hist_string=""):
        if len(correction_df)!=1536:
            print("Applied correction_df is not 1536 long!")
        applied_histo = ROOT.TH2D(histo.GetName()+hist_string, histo.GetName()+hist_string, histo.GetNbinsX(), histo.GetXaxis().GetXmin(), histo.GetXaxis().GetXmax(), histo.GetNbinsY(), histo.GetYaxis().GetXmin(), histo.GetYaxis().GetXmax())
        for n in range(1, histo.GetNbinsX()+1):
            for m in range(1, histo.GetNbinsY()+1):               
                applied_histo.SetBinContent(n, m+int(correction_df.iloc[n-1]*120), histo.GetBinContent(n, m))
        return applied_histo
    
    def delays_to_int(self, df):
        temp_rounded_values = []
        for i, value in enumerate(df["idealDelay"]):
            if value<0:
                temp_rounded_values.append(0)
            else:
                temp_rounded_values.append(round(value, 0))
            #temp_rounded_values.append(math.ceil(value))
            #temp_rounded_values.append(self.general.custom_round(value,1.25)) 
               
        df['bunchDelay'] = temp_rounded_values
        return df
    
    def int_optimizer(self):
        all_delays_df_wInt = {}
        all_int_applied_histos = {}
        fit_stdev_values = {} #For storing the individual stdev values
        int_optimized_df = None
        
        for key in self.all_df_float.keys():
            all_delays_df_wInt[key] = self.delays_to_int(self.all_df_float[key])#DI.delays
            all_int_applied_histos[key] = self.applier(all_delays_df_wInt[key]['bunchDelay'], self.histo, hist_string="_intApplied"+str(key).split(".")[-1])
        
        for i, diffs in enumerate(all_int_applied_histos.items()):
            temp_profileX = diffs[1].ProfileX()
            temp_binContent_values = [] #For storing the bin contents for this histo

            for j in range(0, 1536):
                content = temp_profileX.GetBinContent(j)
                temp_binContent_values.append(content)
                
            temp_stdev = statistics.stdev(temp_binContent_values)
                
            if "stDevs" not in list(fit_stdev_values.keys()):
                fit_stdev_values["stDevs"] = [temp_stdev]
            else:
                fit_stdev_values["stDevs"].append(temp_stdev)
                    
        #Now that we have the correct reference_number_offset index, now we can collect the proper delays from the df
        print(fit_stdev_values)
        min_index = np.argmin(fit_stdev_values["stDevs"])
        min_offset_key = list(all_delays_df_wInt.keys())[min_index] 
        print("Minimum stDev Reference Num: " + min_offset_key)    
        int_optimized_df = all_delays_df_wInt[min_offset_key]

        return int_optimized_df, float(min_offset_key)
    
    def df_to_hist(self, correction_df, pad_df, histo_string=""):
        differences = ROOT.TH1D(self.histo.GetName()+histo_string, self.histo.GetName()+histo_string, self.histo.GetNbinsX(), self.histo.GetXaxis().GetXmin(), self.histo.GetXaxis().GetXmax())
        for i, pid in enumerate(pad_df):
            differences.SetBinContent(differences.FindBin(pid), correction_df[i])
        return differences
    
    def calc_gbt_delay(self):
        gbt = int((self.min_reference_point)*120)# + int((math.ceil(self.int_df.loc[:, 'idealDelay'].mean()) - self.int_df.loc[:, 'idealDelay'].mean())*120)
        temp_df = self.int_df
        temp_df["gbtDelay"] = gbt
        return temp_df
    
    def gbt_applier(self, correction_df, histo, hist_string=""):
        if len(correction_df)!=1536:
            print("Applied correction_df is not 1536 long!")
        applied_histo = ROOT.TH2D(histo.GetName()+hist_string, histo.GetName()+hist_string, histo.GetNbinsX(), histo.GetXaxis().GetXmin(), histo.GetXaxis().GetXmax(), histo.GetNbinsY(), histo.GetYaxis().GetXmin(), histo.GetYaxis().GetXmax())
        for n in range(1, histo.GetNbinsX()+1):
            for m in range(1, histo.GetNbinsY()+1):               
                applied_histo.SetBinContent(n, m-int(correction_df.iloc[n-1]*120/120), histo.GetBinContent(n, m))
        return applied_histo
    
    def format_histos(self):
        self.histo.GetXaxis().SetTitle("Expaned Pad ID")
        self.histo.GetYaxis().SetTitle("Time [bx]")
        self.histo.SetTitle(self.histo.GetName())
        
        self.float_applied_histo.GetXaxis().SetTitle("Expaned Pad ID")
        self.float_applied_histo.GetYaxis().SetTitle("Time [bx]")
        self.float_applied_histo.SetTitle(self.float_applied_histo.GetName())
        
        self.int_applied_histo.GetXaxis().SetTitle("Expaned Pad ID")
        self.int_applied_histo.GetYaxis().SetTitle("Time [bx]")
        self.int_applied_histo.SetTitle(self.int_applied_histo.GetName())
        
        self.gbt_applied_histo.GetXaxis().SetTitle("Expaned Pad ID")
        self.gbt_applied_histo.GetYaxis().SetTitle("Time [bx]")
        self.gbt_applied_histo.SetTitle(self.gbt_applied_histo.GetName())
    
    #For reducing the df from being based on padID to being based on group number 
    def df_reducer_group(self, df):
        temp_df = df.copy()
        return temp_df.drop_duplicates(subset=['fed', 'amc', 'oh', 'vfat', 'group'])
    
    #For reducing the df from being based on padID to being based on gbt number
    def df_reducer_gbt(self, df):
        temp_df = df.copy()
        return temp_df.drop_duplicates(subset=['fed', 'amc', 'oh', 'gbt'])
        
#General functions used in the processing
class generalFunctions():
    def __init__(self):
        pass

    #Iterates over all x bins in a 2d histogram, 
    #fits them with a gaussian and outputs a root file containing:
    #   1) a histo with the means as a function of x bin
    #   2) a histo with the sigmas as a function of the x bin
    #   3) canvases with all of the individual fits
    def fit_2d_histogram(self, h2d, output_file=None, init_params=None, param_limits=None, fit_range=None, max_straddle=False, pol0_from_back=False):
        ROOT.gROOT.SetBatch(True)
        if (fit_range is not None) and (max_straddle is not False):
            print("Warning, you have selected two different fit methods. The first one will be chosen.")
            print("Your options are: Whole histo (default), fit_range, or max_straddle")

        outfile = None
        if output_file is not None:
            outfile = ROOT.TFile(output_file, "RECREATE")

        # Check if the input histogram exists
        if not h2d:
            print("Error: Input histogram not found.")
            if output_file is not None:
                outfile.Close()
            return

        fit_means_hist = ROOT.TH1F("fit_means_"+h2d.GetName(), "Fit Means", h2d.GetNbinsX(), h2d.GetXaxis().GetXmin(), h2d.GetXaxis().GetXmax())
        fit_sigmas_hist = ROOT.TH1F("fit_sigmas_"+h2d.GetName(), "Fit Sigmas", h2d.GetNbinsX(), h2d.GetXaxis().GetXmin(), h2d.GetXaxis().GetXmax())
        canvas = ROOT.TCanvas("canvas", "Fits", 800, 600)

        for binx in range(1, h2d.GetNbinsX()+1):
            h1d = h2d.ProjectionY(f"projection_{binx}", binx, binx)
            h1d_maxbin = h1d.GetMaximumBin()
            h1d_maxX = h1d.GetBinCenter(h1d_maxbin)

            if h1d.GetEntries()==0: #Empty bins cause crashes when calling GetParameter(1)
                continue

            fitted_function = ROOT.TF1("fitted_function", "gaus(0)+pol0(3)")

            if init_params is not None: #Making sure init_params is proper
                if not isinstance(init_params, list):
                    raise TypeError("init_params is not a list")
                elif len(init_params) != 4:
                    raise ValueError("init_params is the incorrect length, 4")
                else:
                    fitted_function.SetParameters(init_params[0], init_params[1], init_params[2])

            if pol0_from_back != False:
                fitted_function.SetParameter(0, h1d.GetBinContent(h1d.GetBin(pol0_from_back)))

            if fit_range is not None:
                fit_min, fit_max = fit_range
                fitted_function.SetRange(fit_min, fit_max)
                fit_option = "QR+"
            elif max_straddle:
                fitted_function.SetRange(0.8*h1d_maxX, 1.2*h1d_maxX)
                fit_option = "QR+"
            else:
                fit_option = "Q"

            if param_limits is not None: #Making sure param_limits is proper
                if not isinstance(param_limits, dict):
                    raise TypeError("param_limits should be a dict of the form {0: [param_max0, param_min0], 1: ..., ...}")
                else:
                    for key in param_limits.keys():
                        fitted_function.SetParLimits(int(key), param_limits[key][0], param_limits[key][1])


            h1d.Fit("fitted_function", fit_option)

            if fitted_function.GetParameter(2) < 0:
                print("Bin ", binx, " has a negative sigma. We will take the absolute value!")

            fit_means_hist.SetBinContent(binx, fitted_function.GetParameter(1))
            fit_sigmas_hist.SetBinContent(binx, abs(fitted_function.GetParameter(2)))

            h1d.Draw()
            fitted_function.Draw("same")
            canvas.Update()

            if outfile is not None:
                canvas.Write(f"fit_canvas_bin_{binx}")

        if outfile is not None:
            fit_means_hist.Write("fit_means_hist")
            fit_sigmas_hist.Write("fit_sigmas_hist")
            fit_means_hist.SetDirectory(0)
            fit_sigmas_hist.SetDirectory(0)
            print("Output made!")
            outfile.Close()

        return fit_means_hist, fit_sigmas_hist

    #For taking an input histogram (1d) and a float to output the differences as another histogram of the same x range/bins as the input.
    def compute_difference_histogram(self, input_hist, referenceNum, hist_name_str=""):
        xmin = input_hist.GetXaxis().GetXmin()
        xmax = input_hist.GetXaxis().GetXmax()
        nbins = input_hist.GetNbinsX()
        output_hist = ROOT.TH1D(input_hist.GetName()+"_differences"+hist_name_str, f"Difference from {referenceNum}", nbins, xmin, xmax)

        for i in range(1, nbins + 1):
            bin_content = input_hist.GetBinContent(i)
            difference = referenceNum - bin_content
            output_hist.SetBinContent(i, difference)
        return output_hist

    #For outputting histo contents as a pandas dataframe
    def histogram_to_df(self, hist, x_label, y_label):
        data = []
        for i in range(1, hist.GetNbinsX() + 1):
            bin_center = hist.GetBinCenter(i)
            bin_value = hist.GetBinContent(i)
            data.append([bin_center, bin_value])

        df = pd.DataFrame(data, columns=[x_label, y_label])
        return df



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
            
            
