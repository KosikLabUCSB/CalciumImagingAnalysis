# Calcium Imaging Analysis

## Description
This repository contains MATLAB scripts for analyzing calcium imaging videos in .tif format. The main script, `calciumImagingAnalysis.m`, processes the videos and generates a report file containing analysis results, as well as a .mat file containing spike train and calcium traces. Additionally, there is a configuration file, `config.m`, where the user can specify directory paths and filenames for input and output to batch process videos.

## Usage

### Prerequisites
Ensure you have the following dependencies installed:
- CaImAn (available on the Flatiron Institute GitHub)
- DSP System Toolbox
- NoRMCorre (available on the Flatiron Institute GitHub)
- Image Acquisition Toolbox
- Image Processing Toolbox
- MATLAB Report Generator
- Parallel Computing Toolbox
- Signal Processing Toolbox
- Statistics and Machine Learning Toolbox

### Setup
1. **Download Dependencies:**
   - CaImAn: You may need to download CaImAn from the Flatiron Institute GitHub .
   - NoRMCorre: Download NoRMCorre from the Flatiron Institute GitHub.

2. **Compilation:**
   - After downloading CaImAn, navigate to the MATLAB command window and run the following command to make the file runnable (Note: Xcode Clang may need to be installed for this step):
     ```
     mex -O -largeArrayDims graph_conn_comp_mex.cpp
     ```

3. **Configuration:**
   - Save a copy of `config.m` as `yourname_config.m`, where 'yourname' represents your name or initials. Ensure to include the underscore. This file will store your directory information.

4. **Update Configuration:**
   - Edit the `yourname_config.m` file and update the paths and filenames as per your directory structure and file names.

### Running the Analysis
- Execute the `calciumImagingAnalysis.m` script and provide the name of your config file as input (enclosed in single quotes). Optionally, you can specify whether to save individual figures as .png files by setting the second argument (saveFigures) to '1'. If not specified, figures will not be saved by default.

Example function calls:
- Save figures: `calciumImagingAnalysis('yourname_config.m', 1);`
- Don't save figures: `calciumImagingAnalysis('yourname_config.m');`

## Output
- **Report File:** A PDF report containing analysis summary and figures will be generated in the specified output reports directory.
- **Data File:** A .mat file containing spike train and deconvolved calcium traces will be saved in the output data directory.
- **Figures (Optional):** If specified, individual figures will be saved as .png files in the output figures directory.

## Notes
- The configuration file (`yourname_config.m`) should be tailored to your specific directory structure and file names.
- Ensure all dependencies are correctly installed and accessible in MATLAB's path.

## Creator and Maintainer
- Creator: Ray Gifford (October 2023)
- Maintainer: Ray Gifford (up until February 2024)
