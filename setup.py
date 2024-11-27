#This script is just to make sure you have the proper folders for generating the data histograms and the delays
#Written by ChatGPT
import os

def ensure_folders_exist(base_folder, subfolders):
    """
    Ensures that the base folder and its specified subfolders exist.
    Creates them if they do not already exist.

    :param base_folder: The base directory path as a string.
    :param subfolders: A list of subfolder paths relative to the base folder.
    """
    # Ensure the base folder exists
    if not os.path.exists(base_folder):
        print(f"Creating base folder: {base_folder}")
        os.makedirs(base_folder)
    else:
        print(f"Base folder already exists: {base_folder}")

    # Ensure each subfolder exists
    for subfolder in subfolders:
        folder_path = os.path.join(base_folder, subfolder)
        if not os.path.exists(folder_path):
            print(f"Creating subfolder: {folder_path}")
            os.makedirs(folder_path)
        else:
            print(f"Subfolder already exists: {folder_path}")

ensure_folders_exist("GEM_mcdonalds_data", [])
ensure_folders_exist("GEM_delays", ["delays", "verification_plots", "verification_plots/final", "verification_plots/initial", "verification_plots/intermediate"])

