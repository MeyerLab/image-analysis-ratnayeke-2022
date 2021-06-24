# image-analysis-ratnayeke-2021
Analysis and plotting code used for Ratnayeke et al. 2021. Included are scripts for generating all figures in paper, as well as the timelapse and fixed microscopy image analysis pipeline used to process imaging experiments. Processed data (single-cell measurements) can be downloaded from Dryad (____) for use in plotting scripts. Plots have been pre-generated and included in folders from the executed scripts. 

Image analysis pipeline has primarly been  used on MCF10A cells, and has been used on HeLa and RPE-1 cells. This pipeline has the following functionaliy:
1) Tracking cells from time-lapse microscopy and quantification of fluorescent reporters (nuclear/cytoplasmic signals, foci analysis)
2) Quantitative image-based cytometry (QIBC) analysis of fixed-cell microscopy, with capabilities to handle multi-round imaging of the same sample.
3) Retrospective Time-lapse Synchronized QIBC (RT-QIBC) analysis of fixed-cell microscopy, using prior information of live-cell tracking measurements of matched cells. 
See Ratnayeke et al. 2021 for methodological details.  

Pipeline is based on Cappell et al., 2016, Cell 166, 167-180. June 30, 2016. 

## Installation
Add */image-analysis-ratnayeke-2021/Tracking_code/Tracking_dependencies/* to MATLAB path. 
## Usage
For recreating plots from Ratnayeke et al. 2021, scripts are included in *Figure_generation/*. Please refer to key.xlsx to find the specific script for a given figure panel. To execute script, change dataDir to the location of the downloaded Dryad data files (path should specify location of .mat files) and then execute script.

For general tracking pipeline usage *Tracking_main/* contains main functions for analysis, and example scripts are included for proper usage in  *Tracking_main/Examples*. For consolidating raw processed single-cell data from pipeline, example code is found in *Example_analysis/*.
