% config.m
% --------------------------------------
% Creator: Ray Gifford, October 2023
% Maintainer: Ray Gifford up until February 2024
%
% This file is used to store the user's personalized directory information
% for running calcium image analysis scripts.
%
% Usage
% - (1) Save this file to a new file called yourname_config.m, where 
%   'yourname' is your name or initials. Remember to include the
%   underscore: The .gitignore file for this repo is configured to ignore 
%   any files containing the string '_config' so it will not impact the 
%   main repo.
% - (2) Update the file with the paths and filenames as indicated below.
%
% This file is then input to the analysis functions, at which time the
%   function loads the needed filenames and paths in order to conduct its
%   analysis.

%% Source path

% Full path to input parent folder
CONFIG.srcDir = '~/Desktop/CalciumSignalingData/2023_12_7_PascaControl';


%% Ouput paths

% Full path to output report files
CONFIG.outputReports= "~/Desktop/CalciumSignalingData/2023_12_7_PascaControl/Reports";

% Full path to output data .mat file (don't include filename)
CONFIG.outputData = "~/Desktop/CalciumSignalingData/2023_12_7_PascaControl/Data";

% Full path to output figures
CONFIG.outputFigures = "~/Desktop/CalciumSignalingData/2023_12_7_PascaControl/Figures";


