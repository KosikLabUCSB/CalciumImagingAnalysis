function [CM, Mdl] = crossvalidateMulti(obj, epochedData, varargin)
% -------------------------------------------------------------------------
% CIA = CalImgAnalysis;
% C = CIA.Classification.crossvalidateMulti(epochedData, classifier, oversample, extractFeatures)
% -------------------------------------------------------------------------
% Creator: Ray Gifford (May 2024)
% Maintainer: Ray Gifford up until May 2024
%
% Input (required)
% - epochedData: epoched time windows for each spike instance from each found
%   neuron
%
% Input (optional)
% - classifier: a string specifying classifier type with name-value pairing ('LDA'
%   [linear discriminate analysis, 'SVM' [support vector machine], 'RF'
%   [random forest]) Default: 'LDA'.
% - oversample: a boolean specifying whether or not to oversample classes
%   with fewer observations, to balances classes. Default: False.
% - extractFeatures: a boolean specifying whether or not to extract
%   additional features for classification: mean, variance, peak-to-peak
%   amplitude, RMS [root mean squared]. Default: false.
% - kernel: a string specifying the kernel type for SVM classification.
%   Default: 'rbf'. Possible inputs: 'linear', 'rbf'
% - gridSearch: a boolean specifying whether or not to perform a grid search
%   over hyperparameters in SVM and RF classification. Default: false.
% - maxEvaluations: number of evaluations to conduct in grid search.
%   Default: 50
%
% Outputs:
% - CM: the confusion matrix from classifcaiton.
% - Mdl: the trained model, for predicting classes
%

ip = inputParser;

% set default values
ip.FunctionName = 'crossvalidateMulti';
ip.addRequired('epochedData',@(x) validateattributes(x,{'numeric'}, {'nonempty'}));
ip.addParameter('classifier', 'LDA', @ischar);
ip.addParameter('oversample', false, @(x) validateattributes(x,{'logical'}, {'nonempty'}));
ip.addParameter('extractFeatures', false, @(x) validateattributes(x,{'logical'}, {'nonempty'}));
ip.addParameter('kernel', 'rbf', @ischar);
ip.addParameter('gridSearch', false, @(x) validateattributes(x,{'logical'}, {'nonempty'}));
ip.addParameter('maxEvaluations', 50, @(x) validateattributes(x,{'numeric'}, {'nonempty'}));


parse(ip,epochedData,varargin{:});

classifier = lower(ip.Results.classifier);
kernel = lower(ip.Results.kernel);
maxEvaluations = ip.Results.maxEvaluations;

%calcium_traces = rand(10, 5, 50); % (10 neurons, 5 spikes, 50 frames)
tempEpochs = epochedData;
tempEpochs(29,:,:)= [];
[nNeurons, nSpikes, nFrames] = size(tempEpochs);

% Reshape the data into 2D matrix
X = reshape(permute(tempEpochs, [2, 1, 3]), nSpikes*nNeurons, nFrames);

% Generate labels for each spike instance
Y = repelem(1:nNeurons, nSpikes)';
labels = linspace(1, nNeurons, nNeurons);

% Capture the class distribution using tabulate
classDist = tabulate(Y);
classCounts = classDist(:, 2); % Second column contains the counts

% Find the maximum class count
[maxCount, maxClass] = max(classCounts);

% Calculate class weights
classWeights = maxCount ./ classCounts; % Inverse of class frequency

if ip.Results.oversample
    
    % Oversample minority classes
    X_oversampled = X;
    Y_oversampled = Y;

    for class = unique(Y)'
        if classCounts(class) < maxCount
            % Calculate how many more samples we need
            numToAdd = maxCount - classCounts(class);
            % Randomly sample with replacement
            idx = find(Y == class);
            addIdx = randsample(idx, numToAdd, true);
            % Add these samples to the dataset
            X_oversampled = [X_oversampled; X(addIdx, :)];
            Y_oversampled = [Y_oversampled; Y(addIdx)];
        end
    end

    X = X_oversampled;
    Y = Y_oversampled;

end

if ip.Results.extractFeatures
    % Extract additional features
    X_features = extractFeatures(X);
    X = [X, X_features];

    % Split data into training and testing sets
    cv = cvpartition(Y, 'HoldOut', 0.3);
    X_train = X(training(cv), :);
    Y_train = Y(training(cv));
    X_test = X(test(cv), :);
    Y_test = Y(test(cv));
end
disp(classifier);
switch classifier
        case 'lda'
        
        % Train the LDA model
        Mdl = fitcdiscr(X_train, Y_train);
        
        case 'svm'
          
            if ip.Results.gridSearch
                % Define the SVM template with a different kernel (e.g., RBF)
                t = templateSVM('KernelFunction', kernel, 'Standardize', true);

                % Specify hyperparameter optimization options
                hyperOpts = struct('AcquisitionFunctionName', 'expected-improvement-plus', ...
                                   'MaxObjectiveEvaluations', maxEvaluations);

                % Train the ECOC model with hyperparameter optimization
                Mdl = fitcecoc(X_train, Y_train, 'Learners', t, ...
                               'OptimizeHyperparameters', 'all', ...
                               'HyperparameterOptimizationOptions', hyperOpts);
            else
                % Define the SVM template with a different kernel (e.g., RBF)
                t = templateSVM('KernelFunction', kernel, 'Standardize', true);
                
                % Train the ECOC model
                Mdl = fitcecoc(X_train, Y_train, 'Learners', t);
            end
            
            
            
        %case 'rf'
   
end


% Predict on the test data
Y_pred = predict(Mdl, X_test);

% Obtain the confusion matrix
confMat = confusionmat(Y_test, Y_pred);

% Plot the confusion matrix
figure;
confusionchart(confMat);
title('Confusion Matrix');
CM = confMat;

% Calculate and display accuracy
accuracy = sum(Y_pred == Y_test) / numel(Y_test);
disp(['Accuracy: ', num2str(accuracy)]);

% Additional metrics (optional)
precision = diag(confMat) ./ sum(confMat, 2);
recall = diag(confMat) ./ sum(confMat, 1)';
f1 = 2 * (precision .* recall) ./ (precision + recall);
disp(['Precision: ', num2str(mean(precision, 'omitnan'))]);
disp(['Recall: ', num2str(mean(recall, 'omitnan'))]);
disp(['F1 Score: ', num2str(mean(f1, 'omitnan'))]);

% Compute the RDM (1 - normalized confusion matrix)
normalizedMatrix = confMat ./ sum(confMat, 2);
rdm = 1 - normalizedMatrix;

% Convert the RDM to a distance matrix
distanceMatrix = pdist(rdm, 'euclidean');
squareDistanceMatrix = squareform(distanceMatrix);

% Perform MDS with a different starting point and criterion
opts = statset('MaxIter', 1000, 'Display', 'final');
[Ymat, stress] = mdscale(squareDistanceMatrix, 2, 'Start', 'random', 'Options', opts, 'Criterion', 'stress');

% Check if MDS was successful
if stress > 1e-5
    warning('MDS did not converge well, stress value: %f', stress);
end

% Plot MDS
figure;
scatter(Ymat(:,1), Ymat(:,2), 100, 'filled');
title('MDS Plot');
xlabel('Dimension 1');
ylabel('Dimension 2');
grid on;
hold on;

% Add class labels to the plot
for i = 1:length(Ymat)
    text(Ymat(i, 1), Ymat(i, 2), string(labels(i)), 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'right');
end
hold off;

% Perform hierarchical clustering
Z = linkage(distanceMatrix, 'average');

% Plot dendrogram
figure;
dendrogram(Z);
title('Dendrogram');
xlabel('Sample Index');
ylabel('Distance');
grid on;

end

function features = extractFeatures(traces)
    nSpikes = size(traces, 1);
    features = zeros(nSpikes, 4); % Assuming 4 additional features

    for i = 1:nSpikes
        trace = traces(i, :);
        features(i, 1) = mean(trace);
        features(i, 2) = var(trace);
        features(i, 3) = max(trace) - min(trace); % Peak-to-peak amplitude
        features(i, 4) = rms(trace); % Root mean square
    end
end