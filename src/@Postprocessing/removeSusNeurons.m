function [vettedTrace,vettedSpikeTrain]  = removeSusNeurons(obj, data, method, threshold)
% -------------------------------------------------------------------------
% CIA = CalImgAnalysis;
% C = CIA.Postprocessing.removeSusNeurons(data, varargin)
% -------------------------------------------------------------------------
% Creator: Ray Gifford (May 2024)
% Maintainer: Ray Gifford up until May 2024
%
% data: struct with traces
% method: method for identifying outliers (e.g., 'baseline', 'frequency', 'correlation', 'peak', 'PCA')
% threshold: threshold value for identifying outliers
    
    vettedTrace = data.trace;
    vettedSpikeTrain = data.spikeTrain;
    tempMatrix = data.rawTrace(:,30:end-5);

    switch method
        case 'baseline'
            % Calculate baseline fluctuation (standard deviation)
            baselineFluctuation = std(tempMatrix, 0, 2);
            % Identify neurons with baseline fluctuation above threshold
            badNeurons = find(baselineFluctuation > threshold);
        case 'frequency'
            % Calculate spiking frequency
            spikeFrequency = sum(diff(tempMatrix > threshold), 2);
            % Identify neurons with low spiking frequency
            badNeurons = find(spikeFrequency < threshold);
        case 'correlation'
           % Compute standard deviation to exclude neurons with low baseline
            baselineFluctuation = std(tempMatrix, 0, 2);
            % Find neurons with low baseline fluctuations
            lowBaselineNeurons = find(baselineFluctuation < max(baselineFluctuation)/10); % 
            
            % Exclude zero baseline neurons from the correlation calculation
            validNeurons = setdiff(1:size(tempMatrix, 1), lowBaselineNeurons);
            validMatrix = tempMatrix(validNeurons, :);
            
            % Compute pairwise correlation coefficients for valid neurons
            correlationMatrix = corr(validMatrix');
            % Average correlation coefficient for each neuron
            avgCorrelation = mean(correlationMatrix, 2);
            % Identify neurons with low average correlation
            badNeurons = validNeurons(avgCorrelation < threshold);
            
        case 'peak'
             % Identify peaks in fluorescence traces for each neuron
            numNeurons = size(vettedTrace, 1);
            numPeaks = zeros(numNeurons, 1);
            for i = 1:numNeurons
                [~, locs] = findpeaks(vettedTrace(i, :));
                numPeaks(i) = length(locs);
            end
            % Identify neurons with fewer peaks than the threshold
            badNeurons = find(numPeaks < threshold);
            
        case 'PCA'
            % Perform PCA
            [~, score, ~, ~, explained] = pca(tempMatrix');
            % Compute contribution of each neuron to principal components
            contribution = score.^2 * (explained ./ sum(explained));
            % Identify neurons with low contribution to dominant components
            badNeurons = find(mean(contribution, 2) < threshold);
        otherwise
            error('Invalid method. Choose from ''baseline'', ''frequency'', ''correlation'', ''peak'', or ''PCA''.');
    end

    % Remove bad neurons from fluorescence matrix
    vettedTrace(badNeurons, :) = [];
    vettedSpikeTrain(badNeurons, :) = [];
    
  
end