function arborDetection(configFile, saveFigures)
% arborDetection.m
% -----------------------------------
% Creator: Ray Gifford (February 2024)
% Maintainer: Ray Gifford up until March 2024
%
% This function analyzes calcium imaging videos in .tif format, and
% outputs a report file containing overview of anaysis, images, and .mat
% file containing spike train and deconvolved calcium traces
%
% Input (required)
% - configFile: Name of your config file, as a string (enclosed in single
%   quotes).
%
% Input (optional)
% - saveFigures: Boolean of whether to save out each individual figure as a .png file. 
%   If empty or not entered, the function will default to FALSE. These
%   figures will otherwise be saved as part of a report file.
%
% Outputs:
% - The function outputs a report file containing summary figures, and a
% data struct containing the deconvolved traces of identified neurons, and
% spike train to a into the "outputReports" and "outputData" directories,
% respectively, as specified in configFile.m. If saveFigures is set to
% TRUE, a .png for every analysis step will be saved into the
% "outputFigures" directory as specified in the configFile.m. If set to
% FALSE, no figures will be saved. Default is FALSE.
%
% Example function calls:
% - Save figures: arborDetection('ray_config.m', 1);
% - Don't save figures: arborDetection('ray_config.m');


% adjustable parameters
% rCAMP vs. gCAMP

% start a local cluster for parallel processing
gcp;

% import report dependency, and add packages
import mlreportgen.report.* % import report package API
addpath(genpath('../../ca_source_extraction-master'));  % add packages to matlab path
addpath(genpath('NoRMCorre-master'));
addpath(genpath('ca_source_extraction'));  % add packages to matlab path
addpath(genpath('NoRMCorre'));

% make sure at least one input was provided
assert(nargin >= 1, 'At least one input (config file name as string) is required.')


% if config filename doesn't end in '.m', append it.
if ~strcmp(configFile(end-1:end), '.m')
    disp('The full config file name with .m, was not specified... but I get what you meant'); 
    configFile = [configFile '.m'];
end

% if saveFigures input not entered or empty, set to 0
if nargin < 2 || isempty(saveFigures)
    saveFigures = 0; end

% Run the config script to get the CONFIG struct.
run(configFile)

addpath(CONFIG.srcDir); % add data folder path
%filenames = ["3_PreKet010", "4_PreKet007", "5_PreKet006", "6_PreKet005", "7_PreKet004", "9_PreKet003", "10_PreKet"]; % create list of all filenames to iterate through

filelist = dir(fullfile(CONFIG.srcDir, '**/*.tif*'));
filenames = string({filelist(:).name});

%% parallel process each video in source folder, in loop

q = 4;
    
% create iterative path
filename = filenames(q);
basePath = CONFIG.srcDir;
path = basePath + "/" + filename;
[folder, baseFileName, extension] = fileparts(path);
disp(path);
disp("Starting video" + baseFileName); % current video file being analyzed


% read file and determine dynamic range
Y = read_file(path);
[d1,d2,T] = size(Y);    % dimensions of file
Y = Y - min(Y(:));      % remove negative offset

minY = quantile(Y(1:9e6),0.0005);
maxY = quantile(Y(1:9e6),1-0.0005);

p = 2;

% data pre-processing
[P,M_nr] = preprocess_data(Y,p);

Cn =  correlation_image(M_nr); %reshape(P.sn,d1,d2);  %max(Y,[],3); %std(Y,[],3); % image statistic (only for display purposes)

% padded selecition
pad = 5;

% test canny edge detector
figure();

im =imagesc(Cn,[min(Cn(:)),max(Cn(:))]);
im2 = double(Cn(pad:d1-pad, pad:d2-pad));
BW1 = edge(im2,'Canny', 0.11);

figure();
imshow(BW1);

% test hough transform
[H,theta,rho] = hough(BW1);
P = houghpeaks(H,70,'threshold',ceil(0.08*max(H(:))));
lines = houghlines(BW1,theta,rho,P,'FillGap',3,'MinLength',1);

hold on
    
max_len = 0;

for k = 1:length(lines)
   xy = [lines(k).point1; lines(k).point2];
   plot(xy(:,1),xy(:,2),'LineWidth',1,'Color','black');
end
 
    
end
