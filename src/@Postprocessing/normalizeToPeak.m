function normalizedEpochedData = normalizeToPeak(obj, epochedData)
% -------------------------------------------------------------------------
% CIA = CalImgAnalysis;
% C = CIA.Classification.clusterEpochs(epochedData, numClusters)
% -------------------------------------------------------------------------
% Creator: Ray Gifford (May 2024)
% Maintainer: Ray Gifford up until May 2024
%
% This function normalizes epoched spike instances for each found neuron,
% to the peak amplitude. This may help to better classify spikes by other
% aspects of their shape, which may better present differences in dynamics,
% rise time, or decay time
%
% epochedData: epoched time windows for each spike instance from each found
% neuron
%

    [nNeurons, nSpikeInstances, nTimePoints] = size(epochedData);

    normalizedEpochedData = zeros(size(epochedData));

    for neuron = 1:nNeurons
        for spike = 1:nSpikeInstances
            spikeData = squeeze(epochedData(neuron, spike, :));
            peakAmplitude = max(spikeData);
            if peakAmplitude > 0
                normalizedEpochedData(neuron, spike, :) = spikeData / peakAmplitude;
            else
                normalizedEpochedData(neuron, spike, :) = spikeData; % Avoid division by zero
            end
        end
    end
end