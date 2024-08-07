function clusterEpochs(obj, epochedData, numClusters)
% -------------------------------------------------------------------------
% CIA = CalImgAnalysis;
% C = CIA.Classification.clusterEpochs(epochedData, numClusters)
% -------------------------------------------------------------------------
% Creator: Ray Gifford (May 2024)
% Maintainer: Ray Gifford up until May 2024
%
% epochedData: epoched time windows for each spike instance from each found
%neuron
% numClusters: number of clusters to find in by unsupervised clustering:
% PCA aided kmeans clustering

    

    [nNeurons, nSpikeInstances, nTimePoints] = size(epochedData);
    
    % Flatten the 3D epochedData into a 2D matrix for PCA
    reshapedData = reshape(epochedData, nNeurons * nSpikeInstances, nTimePoints);
    
    % Perform PCA
    [coeff, score, ~, ~, explained] = pca(reshapedData);
    
    % Determine the number of components to retain 99% variance
    cumulativeVariance = cumsum(explained);
    numComponents = find(cumulativeVariance >= 99, 1);
    
    % Project data onto the retained principal components
    reducedData = score(:, 1:numComponents);
    
    % Conduct k-means clustering
    [clusterIdx, clusterCenters] = kmeans(reducedData, numClusters);
    
    % Map cluster indices back to the original 3D data structure
    clusterIdx3D = reshape(clusterIdx, [nNeurons, nSpikeInstances]);

    % Plot the neurons in the PCA space with cluster coloring
    figure;
    hold on;
    colors = lines(numClusters); % Use 'lines' colormap for distinct cluster colors
    
    if numComponents > 1
        % 2D or higher projection
        for cluster = 1:numClusters
            clusterData = reducedData(clusterIdx == cluster, :);
            scatter(clusterData(:, 1), clusterData(:, 2), 10, 'filled', 'MarkerFaceColor', colors(cluster, :));
        end
        scatter(clusterCenters(:, 1), clusterCenters(:, 2), 100, 'k', 'x', 'LineWidth', 2); % Cluster centers
        xlabel(['Principal Component 1 (' num2str(explained(1), '%.2f') '%)']);
        ylabel(['Principal Component 2 (' num2str(explained(2), '%.2f') '%)']);
    else
        % 1D projection
        for cluster = 1:numClusters
            clusterData = reducedData(clusterIdx == cluster, :);
            scatter(clusterData, zeros(size(clusterData)), 10, 'filled', 'MarkerFaceColor', colors(cluster, :));
        end
        scatter(clusterCenters, zeros(size(clusterCenters)), 100, 'k', 'x', 'LineWidth', 2); % Cluster centers
        xlabel(['Principal Component 1 (' num2str(explained(1), '%.2f') '%)']);
    end
    
    title('Neurons in Proximity Space with K-means Clustering');
    legend(arrayfun(@(x) sprintf('Cluster %d', x), 1:numClusters, 'UniformOutput', false));
    hold off;

    % Plot mean waveforms for each cluster
    plotClusterMeanWaveforms(epochedData, clusterIdx3D, numClusters, nSpikeInstances, nTimePoints);
end

function plotClusterMeanWaveforms(epochedData, clusterIdx3D, numClusters, nSpikeInstances, nTimePoints)
    [nNeurons, ~, ~] = size(epochedData);

    % Initialize figure
    figure;

    % Loop through each cluster
    for cluster = 1:numClusters
        % Get the indices of the neurons in the current cluster
        [neuronIdx, spikeIdx] = find(clusterIdx3D == cluster);

        if isempty(neuronIdx)
            continue;
        end

        % Extract the data for the current cluster
        clusterData = [];
        for i = 1:length(neuronIdx)
            clusterData = [clusterData; squeeze(epochedData(neuronIdx(i), spikeIdx(i), :))'];
        end

        % Calculate the mean waveform for the current cluster
        meanWaveform = mean(clusterData, 1);

        % Plot the mean waveform
        subplot(numClusters, 1, cluster);
        plot(meanWaveform, 'LineWidth', 2);
        title(['Mean Waveform for Cluster ', num2str(cluster)]);
        xlabel('Time Frames');
        ylabel('Amplitude');
    end
end



