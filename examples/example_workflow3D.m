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

filename = 'postKCL_exampleVolume.tif';
Y = read_file(filename);
% Import example data by locating the PGRNKO_Example.tif file on your desktop
%[filename, pathname] = uigetfile(...    
%   {'*.jpg; *.JPG; *.jpeg; *.JPEG; *.img; *.IMG; *.tif; *.TIF; *.tiff, *.TIFF','Supported Files (*.jpg,*.img,*.tiff,)'; ...
%    '*.tif','tif Files (*.tif)';...
%    '*.TIF','TIF Files (*.TIF)';...
%    '*.tiff','tiff Files (*.tiff)';...
%    '*.TIFF','TIFF Files (*.TIFF)'},...    
%    'MultiSelect', 'on');

% Load Data
%path = strcat(pathname,filename);
%[folder, baseFileName, extension] = fileparts(path);
%Y = read_file(path);

% Normalize data
%[P,yNorm] = CIA.Preprocessing.normalize(Y);
%%
% Find all neurons
data = CIA.Analysis.findNeurons(yNorm, P, filename, 'isVolume', true);
%%
% Output data
saveToPath = strcat(folder, '/', baseFileName, '.mat');
save(saveToPath, 'data');

% Plot traces
CIA.Visualization.plotTraces(data.trace, data.spikeTrain);

% Plot spike train
CIA.Visualization.plotSpikeTrain(data.trace, data.spikeTrain, 'heatmap', true);
