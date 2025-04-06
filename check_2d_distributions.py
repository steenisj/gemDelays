#For generating pdfs from the 2d distributions of the timing in the data root files
#Argument 0 is the root_file_path
#Argument 1 is the output_pdf_path
#Argument 2 is optional and gives the canvas string to search for within the root file

import ROOT
from PIL import Image
import os
import sys
import glob
import signal
import shutil

ROOT.gROOT.SetBatch(ROOT.kTRUE)
ROOT.gStyle.SetOptFit(1111)

def hists_2d_to_pdf(directory, output_pdf_path, file_string="GE*_delays.root", hist_string=""):
    # Open the ROOT file
    #root_file = ROOT.TFile.Open(root_file_path)

    # Create a temporary directory to store images
    temp_dir = "./temp_canvas_images"
    if os.path.exists(temp_dir):
       shutil.rmtree(temp_dir)
    os.makedirs(temp_dir, exist_ok=True)

    # List to store image paths
    image_paths = []
    images = []

    #for i, file in enumerate(os.listdir(directory+"/"+file_string)): 
    for i, file in enumerate(glob.glob(directory+"/"+file_string)):
        root_file = ROOT.TFile.Open(file)
        #canvas = ROOT.TCanvas("canvas", "My Canvas", 800, 600)
        #print(root_file)
        keys = root_file.GetListOfKeys()
        if len(keys)==0:
            continue
        key_name_search = file.split("/")[-1].replace("_delays.root","")+hist_string
        for key in keys:
            name = key.GetName()
            if name == key_name_search:
                data_key = name
        
        hist = root_file.Get(data_key)
        hist.RebinY(120)
        print(hist)

        canvas = ROOT.TCanvas("canvas", "My Canvas", 800, 600)
        hist.GetYaxis().SetRangeUser(4,13)
        hist.Draw("COLZ")
        hist.SetTitle(file.split("/")[-1].replace(".root",""))
        image_path = os.path.join(temp_dir, f"fit_2d_hist{i}.png")
        image_paths.append(image_path)
        canvas.SaveAs(image_path)
        root_file.Close()

    for image_path in image_paths:
        images.append(Image.open(image_path))

    if images:
        images[0].save(output_pdf_path, save_all=True, append_images=images[1:])

    print("PDF Output", output_pdf_path)

    # Clean up temporary images
    for image_path in image_paths:
        os.remove(image_path)
    os.rmdir(temp_dir)

def cleanup_and_exit(signum, frame):
    print(f"Received signal {signum}. Cleaning up before exit.")
    os.system('rm -r ./temp_canvas_images')

    # Get the parent process ID (Bash process)
    parent_pid = os.getppid()
    
    # Send SIGTERM to the parent process (Bash script)
    print(f"Killing parent process {parent_pid}")
    os.kill(parent_pid, signal.SIGTERM)

    sys.exit(0)


signal.signal(signal.SIGINT, cleanup_and_exit)
signal.signal(signal.SIGTERM, cleanup_and_exit)


if __name__ == "__main__":
    # Usage
    if len(sys.argv)==3:
        hists_2d_to_pdf(sys.argv[1], sys.argv[2])
    elif len(sys.argv)==5:
        hists_2d_to_pdf(sys.argv[1], sys.argv[2], file_string=sys.argv[3], hist_string=sys.argv[4])

