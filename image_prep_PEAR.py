import pandas as pd
import re
from skimage import io, data, util, restoration
import tifffile as tif
import numpy as np
from skimage.transform import resize
import sys
import glob, os
import csv
import matplotlib.pyplot as plt

# Extract the name of the calcium or voltage image from the parameter file
def extract(x):
    placeholder = str(data[data[0].str.contains(x)])
    placeholder = re.findall('"([^"]*)"', placeholder)
    placeholder = [float(i) for i in placeholder]
    return round(placeholder[0])
def find_calcium_name(x):
    placeholder = data[data[0].str.contains(x)]
    return placeholder

# Select either the Kinetic imaging series or Static image as the output for Cellpose to use for segementation.
#CellposeSelect = "Kinetic"
CellposeSelect = "Static"

# create all of the directories to store the various data in
path = sys.argv[1]
os.chdir(path)
os.makedirs("for_segmentation/", exist_ok=True)
os.makedirs("for_segmentation/std/", exist_ok=True)
os.makedirs("for_segmentation/downsample/", exist_ok=True)
os.makedirs("for_matlab/", exist_ok=True)
for file in glob.glob("*.txt"):
    data = pd.read_csv(file, header = None)


x = find_calcium_name('Acquisition Type = "Kinetic Series"')
calcium_name_full = str(data[0][x.index[-1]-1])
calcium_name = re.findall('"([^"]*)"', calcium_name_full)[0]

# calculate the parameters needed for analysis in MATLAB
fps = extract('Frame Rate')
pace_frequency = extract('Stimulation Pulse Frequency')
delay = extract('Trigger Delay')/1000.0
pace_number = int(extract('Stimulation Pulse Number')-delay*pace_frequency)
total_duration = extract('Series Duration')*(fps-1)
pixels = extract('Height')
pace_duration = int(round(fps*pace_number/pace_frequency))

# write the parameters to a file for later use
parameters = [fps, pace_number, pace_frequency, pixels, total_duration, pace_duration]
with open("for_matlab/parameters.csv", 'w', newline='') as myfile:
     wr = csv.writer(myfile, quoting=csv.QUOTE_NONE)
     wr.writerow(parameters)

os.chdir("scan/")
#****************Downsample kinetic series files and calculate the standard deviation.*****************
for folder in glob.glob("Well_*"):
    print(folder + "/" + str(calcium_name) + '/*.tif')
    im_collection = io.imread_collection(folder + "/" + str(calcium_name) + '/*.tif')
    im_rescaled = np.zeros((len(im_collection),256,256))
    for i in range(0,len(im_collection)):
        im_rescaled[i,:,:] = resize(im_collection[i],(256,256), anti_aliasing=True).astype('float32')
    im_std = np.std(im_rescaled[:pace_duration-1,:,:], axis=0).astype('float32')
    im_std = im_std[np.newaxis,:,:]
    img_min = im_std.min()
    img_max = im_std.max()
    im_std -= img_min
    if img_max > img_min + 1e-3:
        im_std /= (img_max - img_min)
    im_std *= 65536
    tif.imwrite('../for_segmentation/downsample/' + folder + '_downsample.tif', im_rescaled)

# **********Downsample the static image if available.*****************
    print(folder + "/" + "PEAR" + '/*.tif')
    im_collection = io.imread(folder + "/" + "PEAR" + '/*.tif')
    im_rescaled_PEAR = np.zeros((256, 256))
    im_rescaled_PEAR[:, :] = resize(im_collection, (256, 256), anti_aliasing=True).astype('float32')
    im_std2 = im_rescaled_PEAR
    im_std2 = im_std2[np.newaxis, :, :]
    img_min2 = im_std2.min()
    img_max2 = im_std2.max()
    im_std2 -= img_min2
    if img_max2 > img_min2 + 1e-3:
        im_std2 /= (img_max2 - img_min2)

    im_std2 *= 65536
    im_std2 = np.round(im_std2)
    tif.imwrite('../for_segmentation/downsample/' + folder + '_downsample_PEAR.tif', im_rescaled_PEAR)
    background = restoration.rolling_ball(im_std2)
    im_std2 = im_std2 - background

    if CellposeSelect == "Kinetic":
        tif.imwrite('../for_segmentation/std/' + folder + '_std.tif', im_std)
    elif CellposeSelect == "Static":
        tif.imwrite('../for_segmentation/std/' + folder + '_std.tif', im_std2)


