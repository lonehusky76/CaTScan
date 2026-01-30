import pandas as pd
import re
from skimage import io, data, color
import matplotlib.pyplot as plt
import tifffile as tif
import numpy as np
from skimage.transform import rescale, resize, downscale_local_mean
import sys
import glob, os
import csv
import cv2
from skimage.util import img_as_float, img_as_ubyte

def extract(x):
    placeholder = str(data[data[0].str.contains(x)])
    placeholder = re.findall('"([^"]*)"', placeholder)
    placeholder = [float(i) for i in placeholder]
    return round(placeholder[0])
def find_calcium_name(x):
    placeholder = data[data[0].str.contains(x)]
    return placeholder

path = sys.argv[1]
os.chdir(path)
os.mkdir("ap_for_segmentation/")
os.mkdir("ap_for_segmentation/std/")
os.mkdir("ap_for_segmentation/downsample")
for file in glob.glob("*.txt"):
    data = pd.read_csv(file, header = None)

x = find_calcium_name('Acquisition Type = "Kinetic Series"')
calcium_name_full = str(data[0][x.index[0]-1])
calcium_name = re.findall('"([^"]*)"', calcium_name_full)[0]

fps = extract('Frame Rate')
pace_frequency = extract('Stimulation Pulse Frequency')
delay = extract('Trigger Delay')/1000.0
pace_number = int(extract('Stimulation Pulse Number')-delay*pace_frequency)
total_duration = extract('Series Duration')*(fps-1)
pixels = extract('Height')

pace_duration = int(round(fps*pace_number/pace_frequency))

parameters = [fps, pace_number, pace_frequency, pixels, total_duration, pace_duration]
with open("ap_for_matlab/parameters.csv", 'w', newline='') as myfile:
     wr = csv.writer(myfile, quoting=csv.QUOTE_NONE)
     wr.writerow(parameters)


os.chdir("scan/")
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
    im_std *= 255 
    tif.imwrite('../ap_for_segmentation/downsample/' + folder + '_downsample.tif', im_rescaled)
    tif.imwrite('../ap_for_segmentation/std/' + folder + '_std.tif', im_std)




