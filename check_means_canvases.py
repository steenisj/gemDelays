#For generating pdfs from the fit MEANS canvases in the output root files
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

def means_hists_to_pdf(directory, output_pdf_path, file_string="finalFitInformation*.root", hist_string="fit_means_hist", range=None):
    if not range is None:
        joined = ''.join(range)
        cleaned = joined.strip('[]')
        split_values = cleaned.split(',')
        range=split_values

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
        
        hist = root_file.Get(hist_string)
        print(hist)
        canvas = ROOT.TCanvas("canvas", "My Canvas", 800, 600)
        if range is None:
            hist.GetYaxis().SetRangeUser(4,13)
        else:
            hist.GetYaxis().SetRangeUser(float(range[0]),float(range[1]))
        hist.Draw()
        hist.SetTitle(file.split("/")[-1].replace(".root",""))
        image_path = os.path.join(temp_dir, f"{hist_string}{i}.png")
        image_paths.append(image_path)
        canvas.SaveAs(image_path)
        root_file.Close()

    for image_path in image_paths:
        images.append(Image.open(image_path))

    if images:
        images[0].save(output_pdf_path, save_all=True, append_images=images[1:])

    print("PDF Output: ", output_pdf_path)

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
        means_hists_to_pdf(sys.argv[1], sys.argv[2])
    elif len(sys.argv)==4:
        means_hists_to_pdf(sys.argv[1], sys.argv[2], file_string=sys.argv[3])   
    elif len(sys.argv)==6:
        means_hists_to_pdf(sys.argv[1], sys.argv[2], file_string=sys.argv[3], hist_string=sys.argv[4], range=sys.argv[5])
