function epochedData = epochSpikes(obj, fluorescenceMatrix, spikeTrain, minPeakProminence, windowSize)
% -------------------------------------------------------------------------
% CIA = CalImgAnalysis;
% C = CIA.Postprocessing.epochSpikes(fluorescenceMatrix, spikeTrain, minPeakProminence, windowSize)
% -------------------------------------------------------------------------
% Creator: Ray Gifford (May 2024)
% Maintainer: Ray Gifford up until May 2024
%
% flourescenceMatrix: traces
% spikeTrain: spike train
% minPeakProminence: to threshold the peaks, default to include all peaks
% windowSize: size of epoch time windows, default to automatic window
% determination (its not perfect)

    fluorescenceMatrix = fluorescenceMatrix(:,20:end-5);
    spikeTrain = spikeTrain(:,20:end-5);
    [nNeurons, nFrames] = size(fluorescenceMatrix);
    
    % Calculate RMS of the fluorescence data for each neuron
    rmsValues = sqrt(mean(fluorescenceMatrix .^ 2, 2));
    
    % Set default minPeakProminence if not provided
    if nargin < 4 || isempty(minPeakProminence)
        minPeakProminence = 0.01 * rmsValues;
    elseif isscalar(minPeakProminence)
        minPeakProminence = minPeakProminence * ones(nNeurons, 1);  % Apply the same prominence to all neurons
    end
    
    % Determine window size automatically if not provided
    if nargin < 5 || isempty(windowSize)
        % Measure typical spike duration by analyzing the data
        spikeDurations = [];
        for neuron = 1:nNeurons
            % Find spike instances using spikeTrain
            spikeInstances = find(spikeTrain(neuron, :));
            for spike = 1:length(spikeInstances)
                spikeIndex = spikeInstances(spike);
                % Measure spike duration
                spikeEnd = find(fluorescenceMatrix(neuron, spikeIndex:end) < (fluorescenceMatrix(neuron, spikeIndex) / 2), 1, 'first');
                if isempty(spikeEnd)
                    spikeEnd = nFrames - spikeIndex + 1;
                end
                spikeDurations = [spikeDurations, spikeEnd];
            end
        end
        % Set window size to cover typical spike duration
        if ~isempty(spikeDurations)
            avgSpikeDuration = mean(spikeDurations);
            windowSize = ceil(avgSpikeDuration)+5;  % Window size to capture entire spike duration
        else
            windowSize = 10;  % Default window size if no spikes detected
        end
    end
    
    % Preallocate cell array to hold epoched data
    epochedDataCell = cell(nNeurons, 1);
    
    for neuron = 1:nNeurons
        % Find spike instances using spikeTrain
        spikeInstances = find(spikeTrain(neuron, :));
        
        % Skip if no spikes are found for this neuron
        if isempty(spikeInstances)
            continue;
        end
        
        % Preallocate matrix for epoched data for this neuron
        neuronEpochs = zeros(length(spikeInstances), 3 + windowSize);
        
        % Epoch spikes
        for spike = 1:length(spikeInstances)
            spikeIndex = spikeInstances(spike);
            % Ensure window starts 3 frames before the spike onset
            windowStart = spikeIndex - 3;
            windowEnd = spikeIndex + windowSize - 1; % Adjust windowEnd to match windowSize
            
            % Ensure the window does not go out of bounds
            if windowStart < 1
                epochData = [zeros(1, 1 - windowStart), fluorescenceMatrix(neuron, 1:windowEnd)];
                neuronEpochs(spike, 1:length(epochData)) = epochData; % Assign only the length of epochData
            elseif windowEnd > nFrames
                epochData = [fluorescenceMatrix(neuron, windowStart:nFrames), zeros(1, windowEnd - nFrames)];
                neuronEpochs(spike, 1:length(epochData)) = epochData; % Assign only the length of epochData
            else
                epochData = fluorescenceMatrix(neuron, windowStart:windowEnd);
                neuronEpochs(spike, :) = epochData;
            end
        end
        
        % Remove zero or constant epochs
        nonZeroEpochs = ~all(neuronEpochs == 0, 2) & ~all(neuronEpochs == neuronEpochs(:, 1), 2);
        neuronEpochs = neuronEpochs(nonZeroEpochs, :);
        
        % Store the epoched data for this neuron in the cell array
        epochedDataCell{neuron} = neuronEpochs;
    end
    
    % Determine the maximum number of spikes for the most active neuron
    maxSpikes = max(cellfun(@(x) size(x, 1), epochedDataCell));
    
    % Preallocate the 3D matrix for epochedData
    epochedData = nan(nNeurons, maxSpikes, 3 + windowSize); % Use NaN to handle missing data
    
    for neuron = 1:nNeurons
        neuronEpochs = epochedDataCell{neuron};
        if ~isempty(neuronEpochs)
            epochedData(neuron, 1:size(neuronEpochs, 1), :) = neuronEpochs;
        end
    end
    
    % Replace any remaining NaNs with zeros
    epochedData(isnan(epochedData)) = 0;
    
    % Return epochedData with dimensions: nNeurons x nSpikeInstances x (3 + windowSize)
end