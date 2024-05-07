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
load('postKCL_exampleVolume.mat');

% add custom data
% data = uigetfile;
% load(data);

% Instantiate CalImgAnalysis object
CIA = CalImgAnalysis;

% Plot traces
figure();
set(gcf, 'Position', get(0, 'Screensize'));
subplot(1,3,1)
CIA.Visualization.plotTraces(data.rawTrace, data.spikeTrain, 'detectBursts', false, 'chartTitle', 'Raw Traces');

% LPF
fs = 7;
f2=4/fs/2;
[hb2,ha2]=butter(4, f2,'low');
tmpLPF=filtfilt(hb2,ha2,double(full(data.rawTrace)));
trace = tmpLPF;

subplot(1,3,2)
CIA.Visualization.plotTraces(trace, data.spikeTrain,'detectBursts', false, 'chartTitle', 'Low Pass Filtered Traces');

subplot(1,3,3)
CIA.Visualization.plotTraces(data.trace, data.spikeTrain,'detectBursts', false, 'chartTitle', 'Detrended Traces');

% Plot spike train
figure();
CIA.Visualization.plotSpikeTrain(data.trace, data.spikeTrain);