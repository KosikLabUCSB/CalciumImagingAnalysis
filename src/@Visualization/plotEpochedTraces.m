function fig = plotEpochedTraces(obj, epochedData)
% -------------------------------------------------------------------------
% CIA = CalImgAnalysis;
% C = CIA.Visualization.plotEpochedTraces(epochedData)
% -------------------------------------------------------------------------
% Creator: Ray Gifford (May 2024)
% Maintainer: Ray Gifford up until May 2024
%
% EpochedData: the epoched spike instances from each neuron

nNeurons = size(epochedData, 1);
nSpikes = size(epochedData, 2);
nFrames = size(epochedData, 3);



% Number of objects to plot
num_to_plot = 5;

% Create a figure
figure;

% Loop over the first five objects
for neuron = 1:num_to_plot
    
    subplot(num_to_plot, 1, neuron);
    
    hold on;
    
    % Plot each observation for the current object
    for spike = 1:nSpikes
        plot(1:nFrames, squeeze(epochedData(neuron, spike, :)), 'DisplayName', sprintf('Spike %d', spike));
    end
    hold off;
    
    title(sprintf('Neuron %d', neuron));
    xlabel('Time (frames)');
    ylabel('F');
    
    xlim([0 nFrames+5]);
    legend show;
    
end

% Adjust layout
sgtitle('Traces for the First Five Neurons');



end
