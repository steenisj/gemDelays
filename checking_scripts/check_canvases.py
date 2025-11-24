#For generating pdfs from the fit canvases in the output root files
#Argument 0 is the root_file_path
#Argument 1 is the output_pdf_path
#Argument 2 is optional and gives the canvas string to search for within the root file

import ROOT
from PIL import Image
import os
import shutil
import sys
import signal

ROOT.gROOT.SetBatch(ROOT.kTRUE)
ROOT.gStyle.SetOptFit(1111)

def canvases_to_pdf(root_file_path, output_pdf_path, canvas_string="fit_canvas"):
    # Open the ROOT file
    root_file = ROOT.TFile.Open(root_file_path)

    # Create a temporary directory to store images
    temp_dir = "./temp_canvas_images"
    if os.path.exists(temp_dir):
       shutil.rmtree(temp_dir) 
    os.makedirs(temp_dir, exist_ok=True)

    # List to store image paths
    image_paths = []
    images = []

    # Loop over keys in the ROOT file to find canvases
    for key in root_file.GetListOfKeys():
        if key.GetName().startswith(canvas_string):
            canvas_name = key.GetName()
            canvas = root_file.Get(canvas_name)
            canvas.Draw()
            image_path = os.path.join(temp_dir, f"{canvas_name}.png")
            image_paths.append(image_path)
            canvas.SaveAs(image_path)  # Save canvas as PNG image

    # Close the ROOT file
    root_file.Close()

    # Combine images into a single PDF
    for image_path in image_paths:
        images.append(Image.open(image_path))

    if images:
        images[0].save(output_pdf_path, save_all=True, append_images=images[1:])

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
        canvases_to_pdf(sys.argv[1], sys.argv[2])
    elif len(sys.argv)==4:
        canvases_to_pdf(sys.argv[1], sys.argv[2], canvas_string=sys.argv[3])

    print("PDF output: ", sys.argv[2])
