#For generating pdfs from the fit MEANS canvases in the output root files
#Argument 0 is the root_file_path
#Argument 1 is the output_pdf_path
#Argument 2 is optional and gives the canvas string to search for within the root file

import ROOT
from PIL import Image
import os
import sys
import glob

ROOT.gROOT.SetBatch(ROOT.kTRUE)
ROOT.gStyle.SetOptFit(1111)

def means_hists_to_pdf(directory, output_pdf_path, file_string="finalFitInformation*.root"):
    # Open the ROOT file
    #root_file = ROOT.TFile.Open(root_file_path)

    # Create a temporary directory to store images
    temp_dir = "./temp_canvas_images"
    os.makedirs(temp_dir, exist_ok=True)

    # List to store image paths
    image_paths = []
    images = []

    #for i, file in enumerate(os.listdir(directory+"/"+file_string)): 
    for i, file in enumerate(glob.glob(directory+"/"+file_string)):
        root_file = ROOT.TFile.Open(file)
        #canvas = ROOT.TCanvas("canvas", "My Canvas", 800, 600)
        print(root_file)
        hist = root_file.Get("fit_means_hist")
        print(hist)
        canvas = ROOT.TCanvas("canvas", "My Canvas", 800, 600)
        hist.GetYaxis().SetRangeUser(6,13)
        hist.Draw()
        hist.SetTitle(file.split("/")[-1].replace(".root",""))
        image_path = os.path.join(temp_dir, f"fit_means_hist{i}.png")
        image_paths.append(image_path)
        canvas.SaveAs(image_path)
        root_file.Close()

    for image_path in image_paths:
        images.append(Image.open(image_path))

    if images:
        images[0].save(output_pdf_path, save_all=True, append_images=images[1:])

    # Clean up temporary images
    for image_path in image_paths:
        os.remove(image_path)
    os.rmdir(temp_dir)

# Usage
if len(sys.argv)==3:
    means_hists_to_pdf(sys.argv[1], sys.argv[2])
elif len(sys.argv)==4:
    means_hists_to_pdf(sys.argv[1], sys.argv[2], file_string=sys.argv[3])

