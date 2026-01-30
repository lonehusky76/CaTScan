#!/bin/bash

#SBATCH --partition=pcbg-compute               # Uses cardiology node
#SBATCH --account=pcbg
#SBATCH --qos=unlimited
#SBATCH --time=24:00:00                         # Running time (in hours-minute-seconds)  
#SBATCH --mail-type=BEGIN,END,FAIL              # send and email when the job begins, ends or fails
#SBATCH --output=outputs/output_%j.txt                  # Name of the output file
#SBATCH --nodes=1                               # Number of compute nodes
#SBATCH --ntasks=24                             # Number of cpu cores on one node
#SBATCH	--mem=64G					# mem-per-cpu=16G


module load matlab
export CONDA_ENVS_PATH="$HOME/.conda/envs/"

#activates the cellpose environment with Python 3.8.1
source activate cellpose

mkdir $1/for_matlab/
mkdir $1/matlab_output/

#This tests to see if the downsampling has already occurred
#if [-d $1/for_segmentation/ ]; then
#   echo "Directory exists."
#else
#   python image_prep_v2.py $1
#fi

AP_image_prep_XLV.py $1

python -m cellpose --dir $1/for_segmentation/std --pretrained_model Ca_analysis_JM --diameter 0 --verbose --save_tif --no_npy
mv $1/for_segmentation/std/*masks.tif $1/for_segmentation/downsample

#run the script to calculate the roi for each cell and the corresponding voltage and static image value.
python AP_roi_time_XLV.py $1

mv $1/for_segmentation/downsample/*.csv $1/for_matlab

for filename in ${1}/for_matlab/*downsample.csv
do matlab -nodisplay -nosplash -nodesktop -r "addpath('Ca_analysis_matlab_JM'); AP_run_analysis_XLV('${filename}', ['$1','/for_matlab/parameters.csv'],['$1','/matlab_output/'],${2:-0.05},${3:-0.1});exit;"
done

cd $1/matlab_output/
awk '(NR == 1) || (FNR > 1)' *.csv > combined_results.csv


#rm -rf $1/scan/ 
