import skimage
from skimage import io, data, color
from skimage.measure import label, regionprops, regionprops_table
import numpy as np
import glob, os, sys
import pandas as pd 

#Get parameters from command line call of script
#The first parameter is the path of the directory where the files to analyze are. The second is whether to invert the tracing if set to 'TRUE'

if len(sys.argv) < 3:
    WaveInsert = "FALSE"
else:
    WaveInsert = sys.argv[2]

#base_path is the pathway fro the where the downsampled files are.
base_path = sys.argv[1]
os.chdir(base_path)
os.makedirs("for_matlab/static_image_data/", exist_ok=True)

def get_intensity_for_timepoint(intensity_image, label_layer):
    stats = regionprops_table(label_layer, intensity_image=intensity_image, properties=['mean_intensity'])
    return stats['mean_intensity']

def get_intensity(intensity_image_stack, label_layer):
    return [get_intensity_for_timepoint(intensity_image, label_layer) 
            for intensity_image in intensity_image_stack]

i = 0
mask_path = sorted(glob.glob('for_segmentation/downsample/*masks.tif'))
for stack_path in sorted(glob.glob('for_segmentation/downsample/*downsample.tif')):
    mask = io.imread(mask_path[i])
    stack = io.imread(stack_path)
    
    #This either inverts the tracing for Cepheid or not for calcium
    if WaveInsert == 'TRUE':
        x =1/(np.asarray(get_intensity(stack, mask)))
    elif WaveInsert == 'FALSE':
        x =(np.asarray(get_intensity(stack, mask)))

    output_name = os.path.splitext(stack_path)[0]
    DF = pd.DataFrame(x)
    DF.to_csv(output_name + '.csv', index=False)
    i += 1 

#code to measure average fluorescence from static image  channel
k = 0

for stack_path in sorted(glob.glob('for_segmentation/downsample/static_images/*_downsample_PEAR.tif')):
    mask = io.imread(mask_path[k])
    stack = io.imread(stack_path)

    y = np.asarray(get_intensity(stack, mask))

    output_name = os.path.basename(stack_path)
    print('second path ' + output_name)
    TotalCellNumber = len(y.T)

    CellName = [""] * TotalCellNumber
    for j in range(TotalCellNumber):
        CellName[j] = "Cell_ID_" + str(j)

    DF = pd.DataFrame(y.T, CellName)
    DF.columns = ['mean intensity']
    DF.index.name = 'Cell_ID'
    output_name_static_image = base_path + '/for_matlab/static_image_data/' + output_name + '_pos' + '.csv'
    DF.to_csv(output_name_static_image, index=True)

    k += 1

