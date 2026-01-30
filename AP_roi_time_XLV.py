import skimage
from skimage import io, data, color
from skimage.measure import label, regionprops, regionprops_table
import numpy as np
import glob, os, sys
import pandas as pd 

path = sys.argv[1]
os.chdir(path)

def get_intensity_for_timepoint(intensity_image, label_layer):
    stats = regionprops_table(label_layer, intensity_image=intensity_image, properties=['mean_intensity'])
    return stats['mean_intensity']

def get_intensity(intensity_image_stack, label_layer):
    return [get_intensity_for_timepoint(intensity_image, label_layer) 
            for intensity_image in intensity_image_stack]

def find_calcium_name(x):
    placeholder = data[data[0].str.contains(x)]
    return placeholder

#identify the name of the Kinetic Series
 x = find_calcium_name('Acquisition Type = "Kinetic Series"')
        calcium_name_full = str(data[0][x.index[-1] - 1])
        calcium_name = re.findall('"([^"]*)"', calcium_name_full)[0]

# identify the name of the static expression image
        x = find_calcium_name('Acquisition Type = "Single Image"')
        nuc_image_name_full = str(data[0][x.index[-1] - 1])
        nuc_name = re.findall('"([^"]*)"', nuc_image_name_full)[0]

i = 0
mask_path = sorted(glob.glob('for_segmentation/downsample/*masks.tif'))
for stack_path in sorted(glob.glob('for_segmentation/downsample/*downsample.tif')):
    mask = io.imread(mask_path[i])
    stack = io.imread(stack_path)
    x = np.asarray(get_intensity(stack, mask))
    output_name = os.path.splitext(stack_path)[0]
    DF = pd.DataFrame(x)
    DF.to_csv(output_name + '.csv', index=False)
    i += 1 

#code to measure average fluorescence from static image  channel
i = 0
for stack_path in sorted(glob.glob('for_segmentation/downsample/aux/*'+ str(nuc_name) + '.tif')):
    mask = io.imread(mask_path[i])
    stack = io.imread(stack_path)
    x = np.asarray(get_intensity(stack, mask))

    output_name = os.path.splitext(stack_path)[0]
    DF = pd.DataFrame(x)
    DFT = DF.T

    DFT['id'] = DFT.index
    DFT['id'] = 'cell_' + DFT['id'].astype(str)
    DFT = DFT.set_axis(['Mean Fluorescence','Cell ID'], axis=1)
    DFT = DFT.iloc[:,[1,0]]
    DFT.to_csv(output_name + '_pos' + '.csv', index=False)
    i +=1
