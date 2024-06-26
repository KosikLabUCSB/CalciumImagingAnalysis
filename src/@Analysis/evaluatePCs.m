function evaluatePCs(obj, epochedData)
% -------------------------------------------------------------------------
% CIA = CalImgAnalysis;
% C = CIA.Analysis.evaluatePCs(epochedData)
% -------------------------------------------------------------------------
% Creator: Ray Gifford (May 2024)
% Maintainer: Ray Gifford up until May 2024
%
% epochedData: epoched time windows for each spike instance from each found
%neuron


    [nNeurons, nSpikeInstances, nTimePoints] = size(epochedData);

    % Replace NaNs with zeros
    epochedData(isnan(epochedData)) = 0;

    % Flatten the 3D epochedData into a 2D matrix for PCA
    reshapedData = reshape(epochedData, nNeurons * nSpikeInstances, nTimePoints);

    % Perform PCA
    [coeff, score, ~, ~, explained] = pca(reshapedData);
    
    % Get the first three principal components
    numComponents = 3;
    principalComponents = coeff(:, 1:numComponents);
    
    % Plot the loadings for each of the first three principal components
    figure;
    for i = 1:numComponents
        subplot(4, 3, i);
        plot(principalComponents(:, i), 'LineWidth', 2);
        title(['Principal Component ', num2str(i)]);
        xlabel('Time (frames)');
        ylabel('Loading');
    end
    
    % Calculate the mean spike waveform from the filtered data
    meanSpike = mean(reshapedData, 1);

    % Reconstruct synthetic spikes using PC1, PC2, and PC3
    syntheticSpikePC1 = reconstruct_spike(meanSpike, coeff, score, 1);
    syntheticSpikePC2 = reconstruct_spike(meanSpike, coeff, score, 2);
    syntheticSpikePC3 = reconstruct_spike(meanSpike, coeff, score, 3);

    % Plot the original mean spike and synthetic spikes
    subplot(4, 3, 4:6);
    plot(meanSpike, 'b', 'LineWidth', 2);
    hold on;
    if ~isempty(syntheticSpikePC1)
        plot(syntheticSpikePC1, 'r--', 'LineWidth', 2);
    end
    title('Synthetic Spike using PC1');
    xlabel('Time Frames');
    ylabel('Calcium Signal');
    legend('Mean Spike', 'Synthetic Spike (PC1)');
    
    
    subplot(4, 3, 7:9);
    plot(meanSpike, 'b', 'LineWidth', 2);
    hold on;
    if ~isempty(syntheticSpikePC2)
        plot(syntheticSpikePC2, 'r--', 'LineWidth', 2);
    end
    title('Synthetic Spike using PC2');
    xlabel('Time Frames');
    ylabel('Calcium Signal');
    legend('Mean Spike', 'Synthetic Spike (PC2)');
    
    
    subplot(4, 3, 10:12);
    plot(meanSpike, 'b', 'LineWidth', 2);
    hold on;
    if ~isempty(syntheticSpikePC3)
        plot(syntheticSpikePC3, 'r--', 'LineWidth', 2);
    end
    title('Synthetic Spike using PC3');
    xlabel('Time Frames');
    ylabel('Calcium Signal');
    legend('Mean Spike', 'Synthetic Spike (PC3)');
    
    hold off;
    
% Assuming the variable of interest is called 'data' in the .mat file.
data_matrix = cleanedEpochs; % Replace 'cleanedEpochs' with the actual variable name if different

% Reshape the data to 2D
reshaped_data = reshape(data_matrix, size(data_matrix, 1), []);

% Perform PCA
[coeff, score, ~] = pca(reshaped_data);

% Use the first three principal components
pca_data_3d = score(:, 1:3);

% Perform k-means clustering
k = 3; % Assuming 3 clusters
clusters_3d = kmeans(pca_data_3d, k);

% Plot the clustered data in PCA space (3D plot)
figure;
scatter3(pca_data_3d(:, 1), pca_data_3d(:, 2), pca_data_3d(:, 3), 50, clusters_3d, 'filled');
title('K-Means Clustering in 3D PCA Space');
xlabel('PCA1');
ylabel('PCA2');
zlabel('PCA3');
colorbar;
end

function syntheticSpike = reconstruct_spike(meanSpike, coeff, score, pcIndex)
    % Calculate the standard deviation of the scores for the selected PC
    pcScoreStd = std(score(:, pcIndex));
    % Reconstruct the synthetic spike
    syntheticSpike = meanSpike + coeff(:, pcIndex)' * pcScoreStd;
end