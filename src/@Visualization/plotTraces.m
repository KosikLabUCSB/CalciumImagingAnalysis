function fig = plotTraces(obj, trace, spikeTrain, varargin)
% -------------------------------------------------------------------------
% CIA = CalImgAnalysis;
% C = CIA.Visualization.plotTraces(trace, spikeTrain, varargin)
% -------------------------------------------------------------------------
% Creator: Ray Gifford (April 2024)
% Maintainer: Ray Gifford up until May 2024
%
% Given a matrix of calcium traces over time (frames), this function plots the traces for each neuron.  
% Optional name-value pair arugments are described below.
%
% INPUT ARGS (REQUIRED):
%   trace: A calcium flourescence over time 2D matrix (nNeurons x Time (Frames)).
%   spikeTrain: A 2D matrix with neuron spike times (nNeurons x Time
%   (Frames))
%
% INPUT ARGS (OPTIONAL NAME-VALUE PAIRS)
%   'detectBursts': A boolean. TRUE (default) specifies that bursts will be
%       detected and plotted as a vertical line at the time point of
%       detection. FALSE specifies no plotting or detection of bursts.
%   'burstThreshold': A value between 0.1 and 1, where 0.1 specifies 10% of 
%       neurons must first simultaneously to be considered a burst. 0.3 (default)
%       is considered 30% of active neurons fire together for a given point
%       in time. 1 specifies that 100% of acitve neurons must fire together
%       to be detected as a burst.
%   'lineWidth' - Use this parameter to set the width of the vertical lines 
%       used to represent timepoints of bursting. Default 1.
%   'traceColor' - Use this parameter to set the color of the flourescent traces
%       for each neuron. We can either pass in color 
%       abbreviations, full-length color names, or RGB color triplets. Default 'black'.
%   'lineColor' - Use this parameter to set the color of the vertical lines
%       used to represent timpoints of bursting. We can either pass in color 
%       abbreviations, full-length color names, or RGB color triplets. Default 'black'.
%    'fontSize': A number to specify the font size of labels.
%    'traceScale': A number to specify a scalar by which the traces will be
%       multiplied. If flourescence intensity is low, increasing this
%       parameter may help to visualize the traces. Default 2.
%    'traceWidth' - Use this parameter to set the width of the trace lines 
%       used to represent the calcium flourescence over time for each neuron. 
%       Default 1.
%
% OUTPUT ARGS:
%   'img' - Figure corresponding to output plot.

ip = inputParser;

% set default values
ip.FunctionName = 'plotTraces';
ip.addRequired('trace',@ismatrix);
ip.addRequired('spikeTrain',@ismatrix);
ip.addParameter('detectBursts', true, @(x) validateattributes(x,{'logical'}, {'nonempty'}));
ip.addParameter('burstThreshold', 0.3, @(x) isnumeric(x)); 
ip.addParameter('lineWidth', 1, @(x) assert(isnumeric(x), ...
    'lineWidth must be a numeric value'));
ip.addParameter('fontSize', 20, @(x) isnumeric(x));
ip.addParameter('lineColor', 'black');
ip.addParameter('traceColor', [128 128 128]/255);
ip.addParameter('traceScale', 3, @(x) assert(isnumeric(x), ...
    'traceScale must be a numeric value'));
ip.addParameter('traceWidth', 1, @(x) assert(isnumeric(x), ...
    'lineWidth must be a numeric value'));

% parse user inputs
parse(ip, trace, spikeTrain, varargin{:});

% new figure
figure();

% style figure
img = gca;

set(gca,'FontSize',ip.Results.fontSize,'color',[1 1 1]);
img.TitleFontSizeMultiplier = 1.1; % make title larger

% label figure
title('Flourescent traces for each neuron over time');
xlabel('Time (Frames)');
ylabel('Neuron');

hold on;
 
% burst detection

if (ip.Results.detectBursts)
    bursts = sum(spikeTrain,1)/size(spikeTrain,1);
    bursts = (bursts>ip.Results.burstThreshold);
end

% plotting traces and bursts
nNeurons = size(trace,1);

    for i = 1:nNeurons

        tempTrace = trace*ip.Results.traceScale; % scale the flourescent intensity

        plot(1:size(tempTrace,2),tempTrace(i,:)+i, 'Color', ip.Results.traceColor, 'linewidth', ip.Results.traceWidth); % plot detrended trace
        hold on;
        
        ylim([0 size(tempTrace,1)+2]);
        xlim([10 size(tempTrace,2)]); % this shortens the x-axis slightly
        
        hold on; 
        
        if(ip.Results.detectBursts)
            xline(find(bursts), 'Color', ip.Results.lineColor, 'lineStyle', '--', 'linewidth', ip.Results.lineWidth)  % plot potential burst
            legend('Trace', 'Burst');
        end

    end
    
% output
fig = img;

end