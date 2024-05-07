% example_visualization.m
% --------------------------------------
% Example code for visualizing traces and spiketrains
% 
% This script covers the following steps:
%   - Clearing workspace
%   - Instantiating CalImgAnalysis object
%   - Loading example data
%   - Plot Traces
%   - Plot Spike Train
%
% Ray - April, 2024

% Clear workspace
clear all; close all; clc;

% Load example data
load('PGRNKO_Example.mat');

% add custom data
% data = uigetfile;
% load(data);

% Instantiate CalImgAnalysis object
CIA = CalImgAnalysis;

% Unpack data from struct
trace = data.trace;
spikeTrain = data.spikeTrain;

% Plot traces
CIA.Visualization.plotTraces(trace, spikeTrain);

% Plot spike train
CIA.Visualization.plotSpikeTrain(trace, spikeTrain, 'heatmap', true);