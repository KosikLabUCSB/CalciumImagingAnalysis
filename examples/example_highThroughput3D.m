% example_workFlow3D.m
% --------------------------------------
% Example code for using object oriented toolkit. Run this code after
% downloading the example .tif file from the following google drive link: 
% https://drive.google.com/file/d/1IrU5uYg7ab2DziMevCSNk1QsPj1WhZ-T/view?usp=sharing
% You can compare your results with that or the example_visualization.m script.
%
% This script executes the following:
%   - Clear workspace
%   - Instantiate CalImgAnalysis object
%   - Load example data
%   - Preprocess
%   - Find Neurons
%   - Plot Traces
%   - Plot Spike Train
%
% Ray - April, 2024

% Clear console, figures, and workspace
clear all; close all; clc

% Instantiate CalImgAnalysis object
CIA = CalImgAnalysis;



%filename = 'postKCL_exampleVolume.tif';
%Y = read_file(filename);
% Import example data by locating the PGRNKO_Example.tif file on your desktop
srcDir = uigetdir;
filelist = dir(fullfile(srcDir, '**/*.tif*'));
filenames = string({filelist(:).name});


for q = 1:length(filenames)
    
    % create iterative path
    filename = filenames(q);
    basePath = srcDir;
    path = basePath + "/" + filename;
    [folder, baseFileName, extension] = fileparts(path);
    
    %disp("Starting video" + baseFileName); % current video file being analyzed
    %dataTemp = data;
    %filename = 'postKCL_exampleVolume.tif';
%Y = read_file(filename);
% Import example data by locating the PGRNKO_Example.tif file on your desktop
%[filename, pathname] = uigetfile(...    
%   {'*.jpg; *.JPG; *.jpeg; *.JPEG; *.img; *.IMG; *.tif; *.TIF; *.tiff, *.TIFF','Supported Files (*.jpg,*.img,*.tiff,)'; ...
%    '*.tif','tif Files (*.tif)';...
%    '*.TIF','TIF Files (*.TIF)';...
%    '*.tiff','tiff Files (*.tiff)';...
%    '*.TIFF','TIFF Files (*.TIFF)'},...    
%   'MultiSelect', 'on');

% Load Data

    Y = read_file(path);

    % Normalize data
    [P,yNorm] = CIA.Preprocessing.normalize(Y);
    %%
    % Find all neurons
    data = CIA.Analysis.findNeurons(Y, P, 'isVolume', true, 'nSlices', 11);
    %%
    % Output data
    
    % Full path to output data .mat file
    outData= "Data";

    if ~exist(basePath+"/"+outData,'dir')
           mkdir(basePath+"/"+outData)
    end
    saveToPath = strcat(basePath, '/', outData, '/', baseFileName, '.mat');
    save(saveToPath, 'data');

end
