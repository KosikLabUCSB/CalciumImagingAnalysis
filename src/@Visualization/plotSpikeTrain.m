function fig = plotSpikeTrain(obj, trace, spikeTrain, varargin)
% -------------------------------------------------------------------------
% CIA = CalImgAnalysis;
% C = CIA.Visualization.plotSpikeTrain(trace, spikeTrain, varargin)
% -------------------------------------------------------------------------
% Creator: Ray Gifford (April 2024)
% Maintainer: Ray Gifford up until May 2024
%
% Given a matrix containing the neural spike train over time, this function 
% plots the spike train for each neuron.  
% Optional name-value pair arugments are described below.
%
% INPUT ARGS (REQUIRED):
%   trace: A calcium flourescence over time 2D matrix (nNeurons x Time (Frames)).
%   spikeTrain: A 2D matrix with neuron spike times (nNeurons x Time
%   (Frames))
%
% INPUT ARGS (OPTIONAL NAME-VALUE PAIRS)
%   'detectBursts': A boolean. TRUE specifies that bursts will be
%       detected and plotted as a vertical line at the time point of
%       detection. FALSE (default) specifies no plotting or detection of bursts.
%   'heatmap': A boolean. TRUE (default) specifies that each plotted point
%       of the spike train will be heatmapped to the flourescent intensity
%       at that spike time.FALSE specifies no heatmapping of intensities.
%   'popOverlay': A boolean. TRUE (default) specifies that the population
%       cumulative flourescent intensity will overlay the spike train 
%       FALSE specifies no population activity overlay.
%   'burstThreshold': A value between 0.1 and 1, where 0.1 specifies 10% of 
%       neurons must first simultaneously to be considered a burst. 0.3 (default)
%       is considered 30% of active neurons fire together for a given point
%       in time. 1 specifies that 100% of acitve neurons must fire together
%       to be detected as a burst.
%   'lineWidth' - Use this parameter to set the width of the vertical lines 
%       used to represent timepoints of bursting. Default 1.
%       abbreviations, full-length color names, or RGB color triplets. Default 'black'.
%   'lineColor' - Use this parameter to set the color of the vertical lines
%       used to represent timpoints of bursting. We can either pass in color 
%       abbreviations, full-length color names, or RGB color triplets. Default 'black'.
%    'fontSize': A number to specify the font size of labels.
%    'markerWidth' - Use this parameter to set the width of the spike lines 
%       used to represent the a spike instance for each neuron. 
%       Default 1.
%
% OUTPUT ARGS:
%   'img' - Figure corresponding to output plot.

ip = inputParser;

% set default values
ip.FunctionName = 'plotTraces';
ip.addRequired('trace',@ismatrix);
ip.addRequired('spikeTrain',@ismatrix);
ip.addParameter('detectBursts', false, @(x) validateattributes(x,{'logical'}, {'nonempty'}));
ip.addParameter('heatmap', false, @(x) validateattributes(x,{'logical'}, {'nonempty'}));
ip.addParameter('popOverlay', true, @(x) validateattributes(x,{'logical'}, {'nonempty'}));
ip.addParameter('burstThreshold', 0.3, @(x) isnumeric(x)); 
ip.addParameter('lineWidth', 1, @(x) assert(isnumeric(x), ...
    'lineWidth must be a numeric value'));
ip.addParameter('markerWidth', 2, @(x) assert(isnumeric(x), ...
    'lineWidth must be a numeric value'));
ip.addParameter('fontSize', 20, @(x) isnumeric(x));
ip.addParameter('lineColor', 'black');

% parse user inputs
parse(ip, trace, spikeTrain, varargin{:});

% new figure
figure();

% style figure
img = gca;

set(gca,'FontSize',ip.Results.fontSize,'color',[1 1 1]);
img.TitleFontSizeMultiplier = 1.1; % make title larger

% label figure
title('Spike train for each neuron over time');
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
tempTrace = trace;
tempSpikeTrain = spikeTrain;
zerosIndx = tempSpikeTrain(:,:)==0;

for i = 1:nNeurons
     y = linspace(1, length(tempTrace), length(tempTrace));
     
     whenZero = zerosIndx(i,:);
     y(whenZero) = [];
     
     tempSpikeTrain2 = tempSpikeTrain(i,:);
     tempSpikeTrain2(whenZero) = [];
     
     
     if(ip.Results.heatmap)
        cmap = colormap(flip(gray));
        v = rescale(tempTrace, 1, 256); % Nifty trick!
        numValues = length(tempSpikeTrain2);
        markerColors = zeros(numValues, 3);

        % Now assign marker colors according to the value of the data.
        for k = 1 : numValues
            row = round(v(i,y(k)));
            markerColors(k, :) = cmap(row, :);
        end
        
        %markerColors(whenZero') = [];
        
        scatter(y', tempSpikeTrain2*i, [], markerColors, "|", 'LineWidth', ip.Results.markerWidth);
        hold on;

        a = colorbar('eastoutside');
        ylabel(a,'Rescaled Flourescent Intensity','FontSize',14);

    else
        scatter(y', tempSpikeTrain2*i, [], "|", 'Color','black','LineWidth', ip.Results.markerWidth);
        hold on;
     end

    if(ip.Results.popOverlay)
            plot(1:size(tempTrace,2), sum(tempTrace,1), 'red','LineWidth', 1);
            legend('Spike', 'Population Activity');
    end
    
    ylim([0 size(tempTrace,1)+1]);
    xlim([10 size(tempTrace,2)]); % this shortens the x-axis slightly

    hold on; 

end
    
if(ip.Results.detectBursts)
    xline(find(bursts), 'Color', ip.Results.lineColor, 'lineStyle', '--', 'linewidth', ip.Results.lineWidth)  % plot potential burst
    legend('Trace', 'Burst');
end
    
% output
fig = img;

end