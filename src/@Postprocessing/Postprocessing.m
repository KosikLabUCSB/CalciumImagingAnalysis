classdef Postprocessing
% -------------------------------------------------------------------------
% CIA = CalImgAnalysis;
% CIA.Postprocessing
% -------------------------------------------------------------------------
% Ray - May. 22, 2024
%
% This class contains all the postprocessing functions of
% CalImgAnalysis. Units can be checked for quality and reliability, vetted, 
% and epoched for classification or other pipelines.
%
% To call a postprocessing function, you must first create an instance of
% the CalImgAnalysis object, access the Classification member variable, then
% call function of interest.  Below is an example call to the function
% removeSusNeurons:
%
% CIA = CalImgAnalysis;
% N = CIA.Postprocessing.removeSusNeurons(Y);
%
% Below are the list of classification functions contained in this class:
%   - removeSusNeurons()
%   - epochSpikes()
%   - 
%   

   properties
   end
   methods
   end
end