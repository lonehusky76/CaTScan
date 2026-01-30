# CaTScan
End to End pipeline for analysis of high-speed calcium and membrane voltage fluorescent imaging data.

This pipeline performs automated segmentation, analysis, and parameter extraction from high-speed functional imaging data such as from calcium imaging.
The esstenial steps are as follows
1. Data is acquired as multiple TIFF files in a single directory.
2. The Downsample script is run to downsample the data for easier manipulation and data transfer
3. The data is uploaded to a high-performance computing cluster or equivalent
4. the 'AnalysisScript' is modified to point to the appropriate output directory
5. Both Python 3 and Matlab must be installed on the cluster
6. The AnalysisScript is run with the appropriate parameters

The command to run the script is as follows: 
        sbatch AnalysisScript.sh YOUR-DATA-DIRECTORY XXX YY
        Where XXX is the threshold typically 0.05 to 0.08 and YY is the pacing stringency parameter, typically 0.1-0.3
