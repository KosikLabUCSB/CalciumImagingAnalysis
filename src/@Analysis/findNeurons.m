function data = findNeurons(obj, Y, P, filename, varargin)
% -------------------------------------------------------------------------
% CIA = CalImgAnalysis;
% C = CIA.Analysis.findNeurons(Y,P,filename, varargin)
% -------------------------------------------------------------------------
% Creator: Ray Gifford (October 2023)
% Maintainer: Ray Gifford up until May 2024
%
% This function analyzes calcium imaging videos in .tif format, and
% outputs a report file containing overview of anaysis, images, and .mat
% file containing spike train and deconvolved calcium traces
%
% Input (required)
% - Y: .
% - P: .
% - filename: .
%
% Input (optional)
% - isVolume: Boolean of whether to save out each individual figure as a .png file. 
%   If empty or not entered, the function will default to FALSE. These
%   figures will otherwise be saved as part of a report file.
%
% Outputs:
% - .
%
% Example function calls:
%
%

% start a local cluster for parallel processing
gcp;

% import report dependency, and add packages
import mlreportgen.report.* % import report package API
addpath(genpath('../../ca_source_extraction-master'));  % add packages to matlab path
addpath(genpath('NoRMCorre-master'));
addpath(genpath('ca_source_extraction'));  % add packages to matlab path
addpath(genpath('NoRMCorre'));

% make sure at least one input was provided
assert(nargin >= 3, 'At least three inputs (Y, P, and filename) are required.')

ip = inputParser;

% set default values
ip.FunctionName = 'findNeurons';
ip.addRequired('Y');
ip.addRequired('P',@isstruct);
ip.addRequired('filename',@ischar);
ip.addParameter('isVolume', false, @(x) validateattributes(x,{'logical'}, {'nonempty'}));
ip.addParameter('nSlices',[], @(x) validateattributes(x,{'int'}, {'nonempty'}));

% parse user inputs
parse(ip, Y, P, filename, varargin{:});

if ~ip.Results.isVolume

    
    [d1,d2,T1] = size(Y);    % dimensions of file

    % Set parameters for calcium image analysis
        K = 70;                                           % number of components to be found
        tau = 3;                                          % std of gaussian kernel (half size of neuron) 
        %p = 2;

        options = CNMFSetParms(...   
            'd1',d1,'d2',d2,...                         % dimensionality of the FOV
            'p',P.p,...                                   % order of AR dynamics    
            'gSig',tau,...                              % half size of neuron
            'merge_thr',0.91,...                        % merging threshold  
            'nb',2,...                                  % number of background components    
            'min_SNR',3,...                             % minimum SNR threshold
            'space_thresh',0.5,...                      % space correlation threshold
            'cnn_thr',0.2,...                           % threshold for CNN classifier
            'make_avi', 1, 'name', filename);           % can be used to make avi, doesn't work super well

    
    % fast initialization of spatial components using greedyROI and HALS
    [Ain,Cin,bin,fin,center] = initialize_components(double(Y),K,tau,options,P);  % initialize
        
    % update spatial components
    d = d1*d2; 
    Yr = reshape(Y,d,T1);
    [A,b,Cin] = update_spatial_components(Yr,Cin,fin,[Ain,bin],P,options);

    % update temporal components
    P.p = 0;    % set AR temporarily to zero for speed
    [C,f,P,S,YrA] = update_temporal_components(Yr,A,b,Cin,fin,P,options);

    % classify components
    rval_space = classify_comp_corr(Y,A,C,b,f,options);
    ind_corr = rval_space > options.space_thresh;           % components that pass the correlation test
                                            % this test will keep processes

    % further classification with cnn_classifier -- not optional in this script, but may not run. Its ok at classifying
    try  % matlab 2017b or later is needed
        [ind_cnn,value] = cnn_classifier(A,[d1,d2],'cnn_model',options.cnn_thr);
    catch
        ind_cnn = true(size(A,2),1);                        % components that pass the CNN classifier
    end     

    % event exceptionality
    fitness = compute_event_exceptionality(C+YrA,options.N_samples_exc,options.robust_std);
    ind_exc = (fitness < options.min_fitness);

    % select components --- regardless of whether you have CNN setup, move forward
    keep = (ind_corr | ind_cnn) & ind_exc;
    
    A_keep = A(:,keep);
    C_keep = C(keep,:);
    
    % merge found components
    P.merg_thr = 0.95;
    [Am,Cm,K_m,merged_ROIs,Pm,Sm] = merge_components(YrA,A_keep,b,C_keep,f,P,S,options);

    % update
    Pm.p = P.p;    % restore AR value
    [A2,b2,C2] = update_spatial_components(Yr,Cm,f,[Am,b],Pm,options);
    [C2,f2,P2,S2,YrA2] = update_temporal_components(Yr,A2,b2,C2,f,Pm,options);
       
    % refine estimates excluding rejected components
    Pm.p = P.p;    % restore AR value
    [A2,b2,C2] = update_spatial_components(Yr,Cm,f,[Am,b],Pm,options);
    [C2,f2,P2,S2,YrA2] = update_temporal_components(Yr,A2,b2,C2,f,Pm,options);
    
    [A_or,C_or,S_or,P_or] = order_ROIs(A2,C2,S2,P2); % order components
    K_m = size(C_or,1);
    [C_df,~] = extract_DF_F(Yr,A_or,C_or,P_or,options); % extract DF/F values (optional)
    
    
    % detrend fluorescence and extract DF/F values
    df_percentile = 30;
    window = 1000; 

    F = diag(sum(Am.^2))*(C2 + YrA2);  % fluorescence
    Fd = prctfilt(F,df_percentile,window);                      % detrended fluorescence
    Bc = prctfilt((Am'*b)*f2,30,1000,300,0) + (F-Fd);       % background + baseline for each component
    F_dff = Fd./Bc;

    % deconvolve data
    nNeurons = size(F_dff,1);
    C_dec = zeros(size(F_dff));
    S = zeros(size(F_dff));
    kernels = cell(nNeurons,1);
    min_sp = 6;    % find spikes resulting in transients above min_sp x noise level

    for i = 1:nNeurons
        [C_dec(i,:),S(i,:),kernels{i}] = deconvCa(F_dff(i,:), [], min_sp, true, false, [], 20, [], 0);
    end

    centerM = com(A_or,d1,d2);  
    centerM = round(centerM);
    
    data.trace = C_dec;
    data.rawTrace = F_dff;
    data.spikeTrain = S>0;
    data.centers = centerM;
    

    % output mat file with struct
    % save(filename, 'data');
    
% else if this is a volume
else
    

if ~isfield(ip.Results,'nSlices') || isempty(ip.Results.nSlices); nSlices = input('What is the total number of slices? \n'); else nSlices = ip.Results.nSlices; end

assert(~mod(size(Y,3),nSlices), "The number of slices specified does not seem accurate, total number of frames is not divisible by nSlices.");

nFrames = size(Y,3)/nSlices;

Y = double(reshape(Y, [512,512,12,200])); % reshape Y



% Define dimensions and size of dataset

if ndims(Y) == 4
    [d1,d2,d3,T] = size(Y);                            % dimensions of dataset
else
    disp('data was not shaped correctly');
end

d = d1*d2*d3;                                          % total number of pixels

% Set parameters

K = 250;                                          % number of components to be found
tau = [3,3,1];                                    % std of gaussian kernel (size of neuron) 
p = 0;                                            % order of autoregressive system (p = 0 no dynamics for slow imaging rate)
merge_thr = 0.95;                                 % merging threshold

options = CNMFSetParms(...                      
    'd1',d1,'d2',d2,'d3',d3,...                  % dimensions of datasets
    'search_method','dilate',...                 % search locations when updating spatial components
    'maxIter',15,...                             % number of NMF iterations during initialization
    'temporal_iter',2,...                        % number of block-coordinate descent steps 
    'fudge_factor',0.98,...                      % bias correction for AR coefficients
    'merge_thr',merge_thr,...                    % merging threshold
    'gSig',tau,'nb',1 ...
    );
%     'deconv_method','constrained_foopsi',...     % activity deconvolution method

%reshape(P.sn,d1,d2,d3)

[P,Y] = preprocess_data(Y,p);

dims = size(Y);
Cn = correlation_image_3D(double(Y), [], dims); % for large datasets change with reshape(P.sn,d1,d2,d3), %max(Y,[],3); %std(Y,[],3); % image statistic (only for display purposes)

% fast initialization of spatial components using greedyROI and HALS

[Ain,Cin,bin,fin,center] = initialize_components(Y,K,tau,options,P);  % initialize

% display centers of found components
%plotCenteroverY(Cn, center, [d1,d2,d3]);  % plot found centers against max-projections of background image

% update spatial components
Yr = reshape(Y,d,T);
[A,b,Cin] = update_spatial_components(Yr,Cin,fin,[Ain,bin],P,options);

% update temporal components
P.p = 0;
[C,f,P,S,YrA] = update_temporal_components(Yr,A,b,Cin,fin,P,options);

[Am,Cm,K_m,merged_ROIs,Pm,Sm] = merge_components(Yr,A,b,C,f,P,S,options);

[A2,b2,Cm] = update_spatial_components(Yr,Cm,f,[Am,b],Pm,options);
%

[C2,f2,P2,S2,YrA2] = update_temporal_components(Yr,A2,b2,Cm,f,Pm,options);

centerM = com(A2,d1,d2,d3);

if size(centerM,2) == 2
    centerM(:,3) = 1;
end

centerM = round(centerM);
%plotCenteroverY(Cn, centerM, [d1,d2,d3]);  % plot found centers against max-projections of background image

[C_df,~] = extract_DF_F(Yr,A2,C2,P2,options);

% detrend fluorescence and extract DF/F values
df_percentile = 30;
window = 1000; 

F = diag(sum(A2.^2))*(C2 + YrA2);  % fluorescence
Fd = prctfilt(F,df_percentile,window);                      % detrended fluorescence
Bc = prctfilt((A2'*b)*f2,30,1000,300,0) + (F-Fd);       % background + baseline for each component
F_dff = Fd./Bc;

% deconvolve data
nNeurons = size(F_dff,1);
C_dec = zeros(size(F_dff));
S = zeros(size(F_dff));
kernels = cell(nNeurons,1);
min_sp = 4;    % find spikes resulting in transients above min_sp x noise level

for i = 1:nNeurons
        [C_dec(i,:),S(i,:),kernels{i}] = deconvCa(F_dff(i,:), [], min_sp, true, false, [], 20, [], 0);
end

data.trace = C_dec;
data.rawTrace = F_dff;
data.spikeTrain = S>0;
data.centers = centerM;


end



    

end

