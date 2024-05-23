function cleanedData = removeSusEpochs(obj, epochedData, confidence)
% -------------------------------------------------------------------------
% CIA = CalImgAnalysis;
% C = CIA.Postprocessing.removeSusEpochs(epochedData, confidence)
% -------------------------------------------------------------------------
% Creator: Ray Gifford (May 2024)
% Maintainer: Ray Gifford up until May 2024
%
% epochedData: epoched time windows for each spike instance from each found
%neuron
% confidence: the confidence interval by which you would like to determine
% outliers by chi-squared distribution. Default to 75% confidence

    % Replace any remaining zero-filled epochs with NaNs
    %epochedData(epochedData == 0) = NaN;
    
    if nargin < 3
        confidence = 0.75; % Default confidence level
    end

    [numNeurons, numSpikeInstances, numTimePoints] = size(epochedData);
    cleanedData = cell(numNeurons, 1);

    for neuron = 1:numNeurons
        neuronData = squeeze(epochedData(neuron, :, :));

        % Remove rows with zero values across all time points
        neuronData = neuronData(any(neuronData, 2), :);

        % If neuronData is empty after removing zero rows, skip to the next neuron
        if isempty(neuronData)
            cleanedData{neuron} = [];
            continue;
        end

        % Perform PCA
        [coeff, score, ~, ~, explained] = pca(neuronData);
        numComponents = find(cumsum(explained) >= 95, 1); % Adjust variance threshold to 80%

        % Project data onto the retained principal components
        reducedData = score(:, 1:numComponents);

        % Calculate Mahalanobis distance
        mahalanobisDist = mahal(reducedData, reducedData);

        % Determine outliers based on the chi-squared distribution
        chiSquareThreshold = chi2inv(confidence, numComponents); % Adjust confidence level
        isOutlier = mahalanobisDist > chiSquareThreshold;

        % Check if all instances are marked as outliers
        if all(isOutlier)
            % If all instances are outliers, retain the original data or choose a representative subset
            cleanedData{neuron} = neuronData; % Or you can select a subset, e.g., the first instance: neuronData(1, :)
        else
            % Remove outliers
            cleanedData{neuron} = neuronData(~isOutlier, :);
        end
    end

    % Determine the maximum number of non-outlier spikes for the most active neuron
    maxSpikes = max(cellfun(@(x) size(x, 1), cleanedData));

    % Preallocate the 3D matrix for cleaned data
    cleanedDataMatrix = NaN(numNeurons, maxSpikes, numTimePoints);

    for neuron = 1:numNeurons
        neuronEpochs = cleanedData{neuron};
        if ~isempty(neuronEpochs)
            cleanedDataMatrix(neuron, 1:size(neuronEpochs, 1), :) = neuronEpochs;
        end
    end

    cleanedDataMatrix(isnan(cleanedDataMatrix)) = 0;
    
    % Return the cleaned data matrix
    cleanedData = cleanedDataMatrix;
    
    
 [nNeurons, nSpikeInstances, nTimePoints] = size(cleanedData);
    
    % Initialize cleanedData with NaNs
    cleanedData2 = nan(nNeurons, nSpikeInstances, nTimePoints);
    
    for neuron = 1:nNeurons
        neuronData = squeeze(cleanedData(neuron, :, :));
        
        % Find rows that are zero or constant
        zeroRows = all(neuronData == 0, 2);
        constantRows = all(neuronData == neuronData(:, 1), 2);
        
        % Keep rows that are not zero or constant
        cleanedNeuronData = neuronData(~zeroRows & ~constantRows, :);
        
        % If there are no valid rows left, fill with small random noise
       if isempty(cleanedNeuronData)
            cleanedNeuronData = rand(nSpikeInstances, nTimePoints) * 1e-5; % Add small noise
        end
        
        % Ensure the number of rows is consistent with nSpikeInstances
        numValidSpikes = size(cleanedNeuronData, 1);
        if numValidSpikes < nSpikeInstances
            cleanedNeuronData = [cleanedNeuronData; zeros(nSpikeInstances - numValidSpikes, nTimePoints)];
        end
        
        % Assign cleaned data back to the cleanedData matrix
        cleanedData2(neuron, :, :) = cleanedNeuronData;
    end
    
     % Replace any remaining NaNs with small random noise
    cleanedData(isnan(cleanedData)) = rand(sum(isnan(cleanedData(:))), 1) * 1e-5;
    
   
end