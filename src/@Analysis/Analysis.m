classdef Analysis
% -------------------------------------------------------------------------
% CIA = CalImgAnalysis;
% CIA.Analysis
% -------------------------------------------------------------------------
% Ray - April. 28, 2024
%
% This class contains all the classification related functions of
% CalImgAnalysis. Spike trains and traces from classified putative neurons
% generated from classification can be passed to the Visualization module to 
% generate rasters, trace plots, and other figures.
%
% To call a classification function, you must first create an instance of
% the CalImgAnalysis object, access the Classification member variable, then
% call function of interest.  Below is an example call to the function
% findNeurons():
%
% CIA = CalImgAnalysis;
% N = CIA.Analysis.findNeurons(Y);
%
% Below are the list of classification functions contained in this class:
%   - findNeurons()
%   - arborDetection()
%   - 
%   

   properties
   end
   methods
   end
end