# -*- coding: utf-8 -*-
"""
Created on Wed May 28 16:56:19 2025

@author: CH252197
"""


from datetime import datetime


from pathlib import Path


import plotly
import plotly.express as px
import plotly.io as pio
pio.renderers.default = "browser"


import pickle
import os


import NucSeg_SpotDetect
import Config as cf


def imshow3D(data):
    fig = px.imshow(data, animation_frame=0, binary_string=True)
    fig.update_xaxes(showticklabels=False)
    fig.update_yaxes(showticklabels=False)
    fig.update_layout(autosize=False, width=1000, height=1000, coloraxis_showscale=False)
    # Drop animation buttons
    fig['layout'].pop('updatemenus')
    plotly.io.show(fig)
    
    
def file_contains_string(directory, search_string, f_num):
    # directory_path = Path(directory)
    # for file in directory_path.iterdir():
    #     if search_string in file.name:
    #         return True, file.name
    # return False, 'No file'
    directory_path = Path(directory)
    file_list = []
    for file in directory_path.iterdir():
        if search_string in file.name:
            file_list.append( (True, file.name) )
    
    if len(file_list) == 0:    
        return False, 'No file'
    else:
        return file_list[f_num]
    
    
def get_nuc(lab_pc, base_path, force_run, f_num, percentile):
    
    # Determine drive location based on lab_pc flag
    drive_loc = r'G:\Shared drives\Po-Ta_Chen_Projects\Data' if lab_pc else r'H:\Shared drives\Po-Ta_Chen_Projects\Data'
    base_path = drive_loc + base_path  # Construct full path
    
    # Load or Compute Nuclei Data
    # Create a unique identifier for the dataset
    # namestamp = base_path.split('\\')[-3] + '__' + base_path.split('\\')[-1]
    if cf.microscope == 'enders6_sd-confocal':
        namestamp  = base_path.split('\\')[-2]+'__'+base_path.split('\\')[-1]
    else:
        namestamp  = base_path.split('\\')[-3]+'__'+base_path.split('\\')[-1]
    
    
    # Check for precomputed data in a pickle file
    has_file, name = file_contains_string(drive_loc + r'\PoTa-analysis', 'pkl_' + namestamp, f_num)
    
    if has_file & ~force_run:
        # Load existing nuclei data
        with open(drive_loc + r'\PoTa-analysis'+os.sep+name, 'rb') as file:
            nuclei_perFOV = pickle.load(file)
    else:
        
        # Get sorted list of TIFF files from the directory
        folder_path = Path(base_path)
        file_names  = [f.name for f in folder_path.iterdir() if f.is_file()]
        file_names.sort()
        
        # Compute nuclei data if no precomputed file exists
        nuclei_perFOV = NucSeg_SpotDetect.main(base_path, file_names, percentile)
        
        # Save the computed data with a timestamp
        timestamp = datetime.now().strftime('%Y%m%d_%H%M%S')
        with open(drive_loc + fr'\PoTa-analysis\pkl_{namestamp}_{timestamp}.pkl', 'wb') as file:
            pickle.dump(nuclei_perFOV, file)
    
    return nuclei_perFOV, namestamp




