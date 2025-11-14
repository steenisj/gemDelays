#These are the classes that are called when generating the GEM delays.
#Note that all the work is done within the __init__ function of the delayGenerator ... look there first!
#Code by Jacob Steenis, 2024/2025
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
    def __init__(self, histo, histo_name, filename, reference_point=9, rebin_num=8, num_optimize_steps=5):
        if histo==None or histo_name==None:
            self.status = False
        
        else:
            self.status = True
            
            #Getting parameters
            self.general = generalFunctions()
            self.histo = self.hotPadRemover(histo) #Input histo
            self.histo_name = histo_name #Input histo name
            self.filename = filename
            self.reference_point = reference_point
            #self.init_reference_point = init_reference_point
            self.rebin_num = rebin_num
            self.num_optimize_steps = num_optimize_steps

            #To rebin the data into groups!
            self.histo = histo.RebinX(self.rebin_num)
            
            #Extracting the data and generating what the "ideal delays" would be
            print("\033[34m\nGetting data for chamber...\n\033[0m")
            self.station, self.region, self.layer, self.chamber = self.gemPad_stringExtractor()
            self.all_df_float, self.expanded_differences_0, self.means_df, self.all_expanded_difference_hists, self.histo = self.data_generator()            
            self.float_applied_histo = self.applier(self.all_df_float["0.0"]["idealDelay"], self.histo, hist_string="_floatCorrectionApplied")
            
            #Getting the integer delays for the groups
            print("\033[34m\nOptimizing standard deviation from reference point...\n\033[0m")
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
            self.final_amplitudes, self.final_means, self.final_sigmas, self.final_backgrounds = (
                self.general.fit_2d_histogram(self.gbt_applied_histo, 
                    output_file="GEM_delays/verification_plots/final/finalFitInformation_"+self.gbt_applied_histo.GetName()+".root",
                    fit_range=[4,12]
                )
            )
            
            #final = with all information
            self.final_df["padID"] = self.final_df["padID"].astype(int)
            self.final_df["bunchDelay"] = self.final_df["bunchDelay"].astype(int)
            
            #For formatting to be input into the electronics (reducing the final_df)
            self.group_df = self.df_reducer_group(self.final_df)
            self.gbt_df = self.df_reducer_gbt(self.final_df)
    
    def gemPad_stringExtractor(self):
        #pattern = r"gemPad_st(\d+)_R([a-z]+)(\d+)L(\d+)CH(\d+)"
        pattern = r"GE(\d)1_([MP])_(\d+)_L(\d)" #r"GE(\d+)1_([a-z]+)_(\d+)_L(\d+)_data"
        match = re.match(pattern, self.filename.split("/")[-1])
            
        if match:
            station = int(match.group(1))
            region = match.group(2)
            if "M" in region:
                region = int(-1)
            elif "P" in region:
                region = int(1)
            
            layer = int(match.group(4))
            chamber = int(match.group(3))
            return station, region, layer, chamber
        
        else:
            raise ValueError("Incorrect string format for the input!")

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
    
    def data_generator(self):   
        all_expanded_difference_hists = {}
        all_delay_df = {}
        rebinned_hist = self.histo #self.grouper(self.histo)
        fit_amplitudes_hist, fit_means_hist, fit_sigmas_hist, fit_backgrounds_hist = (
            self.general.fit_2d_histogram(
                rebinned_hist, 
                output_file="GEM_delays/verification_plots/initial/fitInformation_"+self.histo_name+".root", 
                fit_range=[2,10]
            )
        )
       
        rebinned_hist = rebinned_hist.Clone()

        # Cycle through a range of reference points between self.reference_point 
        # and self.reference_point+1 to choose whichever removes the bias best        
        for i in np.linspace(0.0, 1.0, self.num_optimize_steps+1)[0:-1]:
            difference_hist = self.general.compute_difference_histogram(
                                        fit_means_hist, 
                                        self.reference_point+i, 
                                        hist_name_str="_"+str(i)
                                    )

            all_expanded_difference_hists[i] = difference_hist
            all_delay_df[str(i)] = self.general.histogram_to_df(difference_hist, "padID", "idealDelay")
            self.apply_info_to_df(all_delay_df[str(i)])          
                    
        return (
            all_delay_df, 
            all_expanded_difference_hists[0], 
            self.general.histogram_to_df(fit_means_hist, "padID", "mean"), 
            all_expanded_difference_hists, 
            rebinned_hist
        )
            
    def applier(self, correction_df, histo, hist_string=""):
        applied_histo = ROOT.TH2D(histo.GetName()+hist_string, 
                                    histo.GetName()+hist_string, 
                                    histo.GetNbinsX(), 
                                    histo.GetXaxis().GetXmin(), 
                                    histo.GetXaxis().GetXmax(), 
                                    histo.GetNbinsY(), 
                                    histo.GetYaxis().GetXmin(), 
                                    histo.GetYaxis().GetXmax()
                                )

        for n in range(1, histo.GetNbinsX()+1):
            for m in range(1, histo.GetNbinsY()+1):              
                #print(n, m)
                applied_histo.SetBinContent(n, m+int(correction_df.iloc[n-1]*120), histo.GetBinContent(n, m))
        return applied_histo
    
    def delays_to_int(self, df):
        temp_rounded_values = []
        for i, value in enumerate(df["idealDelay"]):
            rounded_value = round(value, 0)
            if rounded_value<0:
                temp_rounded_values.append(0)
            else:
                temp_rounded_values.append(rounded_value)
            #temp_rounded_values.append(math.ceil(value))
            #temp_rounded_values.append(self.general.custom_round(value,1.25)) 
               
        df['bunchDelay'] = temp_rounded_values
        return df
    
    def int_optimizer(self):
        all_delays_df_wInt = {}
        all_int_applied_histos = {}
        fit_stdev_values = {} #For storing the individual stdev and mean values
        int_optimized_df = None
        
        for key in self.all_df_float.keys():
            all_delays_df_wInt[key] = self.delays_to_int(self.all_df_float[key])#DI.delays
            all_int_applied_histos[key] = self.applier(all_delays_df_wInt[key]['bunchDelay'], self.histo, hist_string="_intApplied"+str(key).split(".")[-1])

        #Cycle through the different adjustments to the reference number (e.g. +0.2, +0.4, ...)
        for i, corrected_data in enumerate(all_int_applied_histos.items()):
            opt_amplitudes_hist, opt_means_hist, opt_sigmas_hist, opt_backgrounds_hist = (
                self.general.fit_2d_histogram(
                    corrected_data[1], 
                    output_file=f"GEM_delays/verification_plots/intermediate/optimizerCheck_{self.histo_name}_{i}.root", 
                    fit_range=[4,12]
                )
            )

            opt_binContent_values = [] #For storing the bin contents for this histo
            temp_hist = ROOT.TH1D(f"tempHist{i}",f"tempHist{i}",50,7,12)
            for j in range(1, opt_means_hist.GetNbinsX()+1):
                content = opt_means_hist.GetBinContent(j)
                if content and content!=0:
                    temp_hist.Fill(content)
            
            temp_stdev = temp_hist.GetStdDev()
            temp_mean = temp_hist.GetMean()

            if "stDevs" not in list(fit_stdev_values.keys()):
                fit_stdev_values["stDevs"] = [temp_stdev]
            else:
                fit_stdev_values["stDevs"].append(temp_stdev)

            if "means" not in list(fit_stdev_values.keys()):
                fit_stdev_values["means"] = [temp_mean]
            else:
                fit_stdev_values["means"].append(temp_mean)

        # Now that we have the correct reference_number_offset index, now we can collect 
        # the proper delays from the df.
        print("Standard deviation values from optimization: ", fit_stdev_values['stDevs'])
        min_index = np.argmin(fit_stdev_values["stDevs"])
        min_offset_key = list(all_delays_df_wInt.keys())[min_index]
        min_associated_mean = fit_stdev_values["means"][min_index]
        
        print(
                "Minimum standard deviation reference number: ", 
                (self.reference_point + float(min_offset_key)), 
                "bunch crossings \n"
            )

        int_optimized_df = all_delays_df_wInt[min_offset_key]

        #overall_adjustment = float(fit_stdev_values["means"][0]) - self.reference_point

        return int_optimized_df, float(min_offset_key) #, overall_adjustment
    
    def df_to_hist(self, correction_df, pad_df, histo_string=""):
        differences = ROOT.TH1D(self.histo.GetName()+histo_string, 
                                    self.histo.GetName()+histo_string, 
                                    self.histo.GetNbinsX(), 
                                    self.histo.GetXaxis().GetXmin(), 
                                    self.histo.GetXaxis().GetXmax()
                                )

        for i, pid in enumerate(pad_df):
            differences.SetBinContent(differences.FindBin(pid), correction_df[i])
        return differences
    
    def calc_gbt_delay(self):
        gbt = int((self.min_reference_point)*120)
        temp_df = self.int_df
        temp_df["gbtDelay"] = gbt
        return temp_df
    
    def gbt_applier(self, correction_df, histo, hist_string=""):
        applied_histo = ROOT.TH2D(histo.GetName()+hist_string, 
                                    histo.GetName()+hist_string, 
                                    histo.GetNbinsX(), 
                                    histo.GetXaxis().GetXmin(), 
                                    histo.GetXaxis().GetXmax(), 
                                    histo.GetNbinsY(), 
                                    histo.GetYaxis().GetXmin(), 
                                    histo.GetYaxis().GetXmax())

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

    def hotPadRemover(self, hist):
        for bx in range(1, hist.GetNbinsX()+1):
            i_max_content = None
            max_content = 1
            total_X_content = sum(hist.GetBinContent(bx, by) for by in range(1, hist.GetNbinsY() + 1))
           
            if total_X_content < 30:
                continue

            for by in range(1, hist.GetNbinsY()+1):
                if hist.GetBinContent(bx, by) > max_content:
                    max_content = hist.GetBinContent(bx, by)
                    i_max_content = by

            check_chans = [ch for ch in range(1, 16) if abs(ch - int(i_max_content/120)-1) > 2]

            if bx<1400 and bx>1380:
                hist.RebinY(120)
                print(max_content)
                for by in range(1, hist.GetNbinsY()+1):
                    print("BY: ", by, "COUNTS: ", hist.GetBinContent(bx, by))
                
                print(check_chans)
                #break

            for i in check_chans:
                if abs(max_content-hist.GetBinContent(bx, hist.GetYaxis().FindBin(i-1)))/max_content < 0.5:
                    print(f"THROWING AWAY padID {bx} SINCE IT'S HOT!")
                    for by in range(hist.GetNbinsY()):
                        hist.SetBinContent(bx, by, 0)
                        hist.SetBinError(bx, by, 0)
        return hist

#General functions used in the processing
class generalFunctions():
    def __init__(self):
        pass

    # Iterates over all x bins in a 2d histogram, 
    # fits them with a gaussian and outputs a root file containing:
    #   1) a histo with the means as a function of x bin
    #   2) a histo with the sigmas as a function of the x bin
    #   3) canvases with all of the individual fits

    def fit_2d_histogram(self, input_hist, output_file=None, fit_range=None, max_straddle=False):
        h2d = input_hist.Clone()
        #h2d.RebinY(120) #Improves the fitting
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

        fit_amplitudes_hist = ROOT.TH1F("fit_amplitudes_"+h2d.GetName(), 
                                            "Fit Amplitudes", 
                                            h2d.GetNbinsX(), 
                                            h2d.GetXaxis().GetXmin(), 
                                            h2d.GetXaxis().GetXmax()
                                        )

        fit_means_hist = ROOT.TH1F("fit_means_"+h2d.GetName(), 
                                            "Fit Means", 
                                            h2d.GetNbinsX(), 
                                            h2d.GetXaxis().GetXmin(), 
                                            h2d.GetXaxis().GetXmax()
                                    )

        fit_sigmas_hist = ROOT.TH1F("fit_sigmas_"+h2d.GetName(), 
                                            "Fit Sigmas", 
                                            h2d.GetNbinsX(), 
                                            h2d.GetXaxis().GetXmin(), 
                                            h2d.GetXaxis().GetXmax()
                                    )

        fit_backgrounds_hist = ROOT.TH1F("fit_backgrounds_"+h2d.GetName(), 
                                            "Fit Backgrounds", 
                                            h2d.GetNbinsX(), 
                                            h2d.GetXaxis().GetXmin(), 
                                            h2d.GetXaxis().GetXmax()
                                        )

        canvas = ROOT.TCanvas("canvas", "Fits", 800, 600)

        for binx in range(1, h2d.GetNbinsX()+1):
            h1d = h2d.ProjectionY(f"projection_{binx}", binx, binx)
            h1d_maxbin = h1d.GetMaximumBin()
            h1d_maxX = h1d.GetBinCenter(h1d_maxbin)

            if h1d.GetEntries()==0: #Empty bins cause crashes when calling GetParameter(1)
                continue

            fitted_function = ROOT.TF1("fitted_function", "gaus(0)+pol0(3)")

            fitted_function.SetParameters(h1d.GetMaximum(), h1d.GetMean(), 0.5)

            if fit_range is not None:
                fit_min, fit_max = fit_range
                fitted_function.SetRange(fit_min, fit_max)
                fit_option = "QR+"
            elif max_straddle:
                fitted_function.SetRange(0.8*h1d_maxX, 1.2*h1d_maxX)
                fit_option = "QR+"
            else:
                fit_option = "Q"

            fitted_function.SetParLimits(0, h1d.GetMaximum()*0.8, h1d.GetMaximum()*1.2)
            fitted_function.SetParLimits(1, h1d.GetMean()*0.3, h1d.GetMean()*1.7)
            fitted_function.SetParLimits(2, 0.1, 4)
            fitted_function.SetParLimits(3, 0, 100)
            
            h1d.Fit("fitted_function", fit_option)
            if fitted_function.GetParameter(2) < 0:
                print("Bin ", binx, " has a negative sigma. We will take the absolute value!")
            
            fit_amplitudes_hist.SetBinContent(binx, fitted_function.GetParameter(0))
            fit_means_hist.SetBinContent(binx, fitted_function.GetParameter(1))
            fit_sigmas_hist.SetBinContent(binx, abs(fitted_function.GetParameter(2)))
            fit_backgrounds_hist.SetBinContent(binx, fitted_function.GetParameter(3))

            h1d.Draw()
            fitted_function.Draw("same")
            canvas.Update()
            
            if outfile is not None:
                canvas.Write(f"fit_canvas_bin_{binx}")

        if outfile is not None:
            fit_amplitudes_hist.Write("fit_amplitudes_hist")
            fit_means_hist.Write("fit_means_hist")
            fit_sigmas_hist.Write("fit_sigmas_hist")
            fit_backgrounds_hist.Write("fit_backgrounds_hist")
            fit_amplitudes_hist.SetDirectory(0)
            fit_means_hist.SetDirectory(0)
            fit_sigmas_hist.SetDirectory(0)
            fit_backgrounds_hist.SetDirectory(0)
            print(f"Output {outfile.GetName()} made!")
            outfile.Close()

        return fit_amplitudes_hist, fit_means_hist, fit_sigmas_hist, fit_backgrounds_hist

    # For taking an input histogram (1d) and a float to output the differences as 
    # another histogram of the same x range/bins as the input.
    def compute_difference_histogram(self, input_hist, referenceNum, hist_name_str=""):
        xmin = input_hist.GetXaxis().GetXmin()
        xmax = input_hist.GetXaxis().GetXmax()
        nbins = input_hist.GetNbinsX()
        output_hist = ROOT.TH1D(input_hist.GetName()+"_differences"+hist_name_str, 
                                    f"Difference from {referenceNum}", 
                                    nbins, 
                                    xmin, 
                                    xmax
                                )

        prev_non_empty_value = 0 #For fixing empty bins
        for i in range(1, nbins + 1):
            bin_content = input_hist.GetBinContent(i)
            difference = referenceNum - bin_content

            if bin_content != 0:
                output_hist.SetBinContent(i, difference)
                prev_non_empty_value = 0
            else: #For fixing empty bins
                output_hist.SetBinContent(i, prev_non_empty_value)

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
            if "GE" in str(key.GetName()):
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
            
            
