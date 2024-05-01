function calciumImagingAnalysis(configFile, saveFigures)
% calciumImagingAnalysis.m
% -----------------------------------
% Creator: Ray Gifford (October 2023)
% Maintainer: Ray Gifford up until February 2024
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
% - Save figures: calciumImagingAnalysis('ray_config.m', 1);
% - Don't save figures: calciumImagingAnalysis('ray_config.m');


% adjustable parameters


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
    disp('You did not inude the full file name... but I understood what you meant'); 
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


%outputReportPath = CONFIG.outputReports;
%outputDataPath = CONFIG.outputData;
%outputFiguresPath = CONFIG.outputFigures;

%% parallel process each video in source folder, in loop

for q = 1:length(filenames)
    
    % create iterative path
    filename = filenames(q);
    basePath = CONFIG.srcDir;
    path = basePath + "/" + filename;
    [folder, baseFileName, extension] = fileparts(path);
    disp(path);
    disp("Starting video" + baseFileName); % current video file being analyzed
    
    % initiate output report file for current video
    rpt = Report(CONFIG.outputReports + "/" + baseFileName,'pdf');
    add(rpt,TitlePage('Title',"Pictorial Overview of Analysis for Video: "+baseFileName,'Author','Ray Gifford'));
    add(rpt,TableOfContents);
    
    % read file and determine dynamic range
    Y = read_file(path);
    [d1,d2,T] = size(Y);    % dimensions of file
    Y = Y - min(Y(:));      % remove negative offset

    minY = quantile(Y(1:9e6),0.0005);
    maxY = quantile(Y(1:9e6),1-0.0005);
   
    % perform motion correction (start with rigid) -- This shouldn't do much if the slice is stuck nicely to the glass
    % parameters motion correction
    % 'd1','d2': size of FOV
    % 'bin_width': how often to update the template
    % 'max_shift': maximum allowed rigid shift

    options_rg = NoRMCorreSetParms('d1',d1,'d2',d2,'bin_width',50,'max_shift',10);

    [M_rg,shifts_rg,template_rg] = normcorre_batch(Y,options_rg);

    
    % downsample video
    tsub = 5;   % downsampling factor (only for display purposes)
    Y_sub = downsample_data(Y,'time',tsub);

    
    % write downsampled video file
     v = VideoWriter("downSampled_video_AVI", 'Grayscale AVI');
     open(v)
     writeVideo(v,mat2gray(Y_sub));
     
     
    % Set parameters for calcium image analysis
        K = 70;                                           % number of components to be found
        tau = 3;                                          % std of gaussian kernel (half size of neuron) 
        p = 2;

        options = CNMFSetParms(...   
            'd1',d1,'d2',d2,...                         % dimensionality of the FOV
            'p',p,...                                   % order of AR dynamics    
            'gSig',tau,...                              % half size of neuron
            'merge_thr',0.91,...                        % merging threshold  
            'nb',2,...                                  % number of background components    
            'min_SNR',3,...                             % minimum SNR threshold
            'space_thresh',0.5,...                      % space correlation threshold
            'cnn_thr',0.2,...                           % threshold for CNN classifier
            'make_avi', 1, 'name', filename);           % can be used to make avi, doesn't work super well
    
    % data pre-processing
    [P,M_nr] = preprocess_data(M_rg,p);
    
    % fast initialization of spatial components using greedyROI and HALS
    [Ain,Cin,bin,fin,center] = initialize_components(double(M_nr),K,tau,options,P);  % initialize

    % output figure displaying centers of found components to report file
    ch = Chapter('ROI Centroids After Fast Initialization'); % create report chapter
    
    Cn =  correlation_image(M_nr); %reshape(P.sn,d1,d2);  %max(Y,[],3); %std(Y,[],3); % image statistic (only for display purposes)
    
    f1 = figure('visible', 'off');
        imagesc(Cn);
        axis equal; axis tight; hold all;
        scatter(center(:,2),center(:,1),'mo');
        title('Center of ROIs found from initialization algorithm');
        drawnow;
    
    fig1 = Figure(f1);
    fig1.Snapshot.Caption = 'Test Caption for this Figure';
    fig1.Snapshot.Height = '7in';
    fig1.Snapshot.Width = '7in';
    fig1.Snapshot.ScaleToFit = true;
    
    add(ch,fig1); % add figure to chapter
    add(rpt,ch);        % add chapter to report
    
    
    % update spatial components
    [d1,d2,T1] = size(M_nr);
    d = d1*d2; 
    Yr = reshape(M_nr,d,T1);
    [A,b,Cin] = update_spatial_components(Yr,Cin,fin,[Ain,bin],P,options);

    % update temporal components
    P.p = 0;    % set AR temporarily to zero for speed
    [C,f,P,S,YrA] = update_temporal_components(Yr,A,b,Cin,fin,P,options);

    % classify components
    rval_space = classify_comp_corr(M_nr,A,C,b,f,options);
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

    % display kept and discarded components
    ch = Chapter('Kept and Discarded Components'); % create report chapter
    
    A_keep = A(:,keep);
    C_keep = C(keep,:);
    
    f2 = figure('visible', 'off');
        subplot(121); montage(extract_patch(A(:,keep),[d1,d2],[30,30]),'DisplayRange',[0,0.15]);
            title('Kept Components');
        subplot(122); montage(extract_patch(A(:,~keep),[d1,d2],[30,30]),'DisplayRange',[0,0.15])
            title('Discarded Components');
    
    fig2 = Figure(f2);
    %fig2.Snapshot.Caption = 'Test Caption for this Figure';
    fig2.Snapshot.Height = '7in';
    fig2.Snapshot.Width = '7in';
    
    add(ch,fig2);       % add figure to chapter
    add(rpt,ch);        % add chapter to report
    
    % merge found components
    P.merg_thr = 0.95;
    [Am,Cm,K_m,merged_ROIs,Pm,Sm] = merge_components(YrA,A_keep,b,C_keep,f,P,S,options);

    % update
    Pm.p = p;    % restore AR value
    [A2,b2,C2] = update_spatial_components(Yr,Cm,f,[Am,b],Pm,options);
    [C2,f2,P2,S2,YrA2] = update_temporal_components(Yr,A2,b2,C2,f,Pm,options);
    
    % display merged component example
    ch = Chapter('Temporal Merging'); % create report chapter
    
    display_merging = 1; % flag for displaying merging example
    
    if and(display_merging, ~isempty(merged_ROIs))
        i = 1; %randi(length(merged_ROIs));
        ln = length(merged_ROIs{i});
        
        f3 = figure('visible', 'off');
            set(gcf,'Position',[300,300,(ln+2)*300,300]);
            for j = 1:ln
                subplot(1,ln+2,j); imagesc(reshape(A_keep(:,merged_ROIs{i}(j)),d1,d2)); 
                    title(sprintf('Component %i',j),'fontsize',16,'fontweight','bold'); axis equal; axis tight;
            end
            subplot(1,ln+2,ln+1); imagesc(reshape(Am(:,K_m-length(merged_ROIs)+i),d1,d2));
                    title('Merged Component','fontsize',16,'fontweight','bold');axis equal; axis tight; 
            subplot(1,ln+2,ln+2);
                plot(1:T1,(diag(max(C_keep(merged_ROIs{i},:),[],2))\C_keep(merged_ROIs{i},:))'); 
                hold all; plot(1:T1,Cm(K_m-length(merged_ROIs)+i,:)/max(Cm(K_m-length(merged_ROIs)+i,:)),'--k')
                title('Temporal Components','fontsize',16,'fontweight','bold')
            drawnow;
            
            fig3 = Figure(f3);
            fig3.Snapshot.Caption = 'Test Caption for this Figure';
            fig3.Snapshot.Height = '7in';
            fig3.Snapshot.Width = '7in';
            add(ch,fig3);       % add figure to chapter
            add(rpt,ch);        % add chapter to report
    end
    
    
    
    % refine estimates excluding rejected components
    Pm.p = p;    % restore AR value
    [A2,b2,C2] = update_spatial_components(Yr,Cm,f,[Am,b],Pm,options);
    [C2,f2,P2,S2,YrA2] = update_temporal_components(Yr,A2,b2,C2,f,Pm,options);
    
    % plot contour
    ch = Chapter('Contour of Spatial Footprints'); % create report chapter
    
    [A_or,C_or,S_or,P_or] = order_ROIs(A2,C2,S2,P2); % order components
    K_m = size(C_or,1);
    [C_df,~] = extract_DF_F(Yr,A_or,C_or,P_or,options); % extract DF/F values (optional)
    
    
    
    f4 = figure('visible', 'off');
    [Coor,json_file] = plot_contours(A_or,Cn,options,1); % contour plot of spatial footprints
    %savejson('jmesh',json_file,'filename');        % optional save json file with component coordinates (requires matlab json library)
    
    fig4 = Figure(f4);
    %fig4.Snapshot.Caption = 'Test Caption for this Figure';
    fig4.Snapshot.Height = '7in';
    fig4.Snapshot.Width = '7in';
    
    add(ch,fig4); % add figure to chapter
    add(rpt,ch);        % add chapter to report
    
    
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
    
    
    % plot a random component
    ch = Chapter('Traces of Example Random Component'); % create report chapter
    i = randi(nNeurons);

    f5 = figure('visible', 'off');
    plot(1:T,F_dff(i,:),'--k'); %plot DF/F trace
    hold all; 
    plot(1:T,C_dec(i,:),'r','linewidth',2); %plot deconvoled signal
    spt = find(S(i,:));

    if spt(1) == 1 
       spt(1) = []; 
    end

    hold on; 
    scatter(spt,repmat(-0.25,1,length(spt)),'m*') % plot spikes
    title(['Component ',num2str(i)]);
    legend('Fluorescence DF/F','Deconvolved','Spikes')
    
    fig5 = Figure(f5);
    %fig5.Snapshot.Caption = 'Test Caption for this Figure';
    fig5.Snapshot.Height = '7in';
    fig5.Snapshot.Width = '7in';
    
    add(ch,fig5); % add figure to chapter
    add(rpt,ch);        % add chapter to report
    
    
    % plot deconvolved traces and spike raster overlayed
    ch = Chapter('Traces and Raster for All Found Components'); % create report chapter
    f6 = figure('visible','off');
    g = gca;

    xlim([100 T]);             % change limit for X to exclude start of recording
     % change limit for Y to trim to number of chosen neurons

    title('Deconvolved Traces and Spike Raster Overlay');
    xlabel('Time (Frames)');
    ylabel('Neuron');
    f6.Position = [200 200 2000 1000];
    %g.FontSize = 20;
    g.TitleFontSizeMultiplier = 1.1;

    hold on;

    for i = 1:nNeurons

        tempDec = C_dec*3;
        tempS = S;
        spt = find(tempS(i,:));


        plot(1:T,tempDec(i,:)+i,'b','linewidth',1); %plot deconvoled signal
        hold on;
        scatter(spt,repmat(-0.25,1,length(spt))+i, 18,'r*') % plot spikes
        hold on;
        xline(115, '--k')  % for potential burst marker
       
        ylim([0 size(tempS,1)+2]);
        hold on;
        legend('Deconvolved','Spike', 'Burst');

    end
    fig6 = Figure(f6);
    fig6.Snapshot.Caption = 'Test Caption for this Figure';
    fig6.Snapshot.Height = '10in';
    fig6.Snapshot.Width = '10in';
    add(ch,fig6); % add figure to chapter
    add(rpt,ch);        % add chapter to report

    % create output struct
    data.spikeTrain = S;
    data.trace = C_dec;

    % output mat file with struct
    outputPath = CONFIG.outputData+ "/" + baseFileName;
    save(outputPath, 'data');
    
    % close report to save
    close(rpt);
end
disp("Finished analysis of all videos in source folder");
