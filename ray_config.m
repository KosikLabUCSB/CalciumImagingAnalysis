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

%% Source path -- Change this to your own source directory (where your videos are)

srcDir = "/MATLAB Drive/2023_11_20_PreKetGreg";

% Full path to input parent folder
CONFIG.srcDir = srcDir;

%% Ouput paths

% Full path to output report files
outReports = "Reports";

if ~exist(srcDir+"/"+outReports,'dir')
       mkdir(srcDir+"/"+outReports)
       CONFIG.outputReports= srcDir + "/" + outReports;

else
    CONFIG.outputReports= srcDir + "/" + outReports;
end
    

% Full path to output data .mat file
outData= "Data";

if ~exist(srcDir+"/"+outData,'dir')
       mkdir(srcDir+"/"+outData)
       CONFIG.outputData= srcDir + "/" + outData;

else
    CONFIG.outputData= srcDir + "/" + outData;
end

% Full path to output figures
outFigs = "Figures";

if ~exist(srcDir+"/"+outFigs,'dir')
       mkdir(srcDir+"/"+outFigs)
       CONFIG.outputFigures= srcDir + "/" + outFigs;

else
    CONFIG.outputFigures= srcDir + "/" + outFigs;
end


