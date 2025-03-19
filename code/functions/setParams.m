% Configures various predefined parameters and saves it into a structure 
% variable 'cfg'. 
function params = setParams(header)
    %%%%%%%%%%%%%%%%%%%%%%%%%
    %% General Information %%
    %%%%%%%%%%%%%%%%%%%%%%%%%
    params.fsamp = header.SampleRate;
    load chanlocs64_2.mat %loads as chanlocs
    chanRemove = {'FP1','FP2','FPZ','M1','M2','EOG'};
    chanIndices = find(ismember(upper({chanlocs.labels}),chanRemove));
    chanlocs(chanIndices) = [];
    [keepIdx, eegIndex] = ismember(header.Label, upper({chanlocs.labels}));
    params.chanlocs = chanlocs(eegIndex(keepIdx));    

    params.eegChannels = 1:64; 
    params.eogChannels = 65:66;
    params.triggerChannel = 67;

    params.chanLabels = header.Label;
    params.chanLabels(65:67)=[];

    params.channelPlot = find(strcmp({params.chanlocs.labels}, 'PO8'));  % normally this one
    params.plotOption = {'LineWidth', 2};
    params.plotColor = {
    [1, 0, 0],        % Red
    [0, 0, 1],        % Blue
    [1, 0.6, 0.6],    % Light Red
    [0.6, 0.6, 1],    % Light Blue
    [0, 1, 0],         % Green
    [0, 0, 0]         % Black
            };
    
    %%%%%%%%%%%%%%
    %% Epoching %%
    %%%%%%%%%%%%%%
    params.epochSamples = -0.5*params.fsamp+1:1.0*params.fsamp;
    params.epochTime = params.epochSamples./params.fsamp;
    params.epochOnset = find(params.epochTime == 0);
    
    %%%%%%%%%%%%%%%%%%%%%
    %% Epoch Rejection %%
    %%%%%%%%%%%%%%%%%%%%%
    params.epochRejection.isCompute = true;
    % params.epochRejection.time = round(0.2*params.fsamp)+1:round(0.8*params.fsamp);
    params.epochRejection.time = round(0.15*params.fsamp)+1:round(0.5*params.fsamp);
    params.epochRejection.time = params.epochRejection.time + params.epochOnset;
    
    %%%%%%%%%%%%%%%%%%%%%
    %% Channel Removal %%
    %%%%%%%%%%%%%%%%%%%%%
    params.channelRemoval.isCompute = false;
    
    %%%%%%%%%%%%%%%%%%%%%
    %% Spectral Filter %%
    %%%%%%%%%%%%%%%%%%%%%
    params.spectralFilter.freqs = [1 30];  % cut-off frequencies
    params.spectralFilter.order = 2;  % 2*params.fsamp for FIR filter
    
    params.EOG.spectralFilter.freqs = [1 10];  % cut-off frequencies
    params.EOG.spectralFilter.order = 2;  % 2*params.fsamp for FIR filter

    %%%%%%%%%%%%%%%%%%%%
    %% Spatial Filter %%
    %%%%%%%%%%%%%%%%%%%%
    params.spatialFilter.type = 'None';  % Option : CAR, Laplace, xDAWN, CCA, CSD, None
    params.spatialFilter.time = round(0.15*params.fsamp)+1:round(0.5*params.fsamp);
    params.spatialFilter.time = params.spatialFilter.time + params.epochOnset;
    params.spatialFilter.nComp = 3;
    params.spatialFilter.classes = [1, 2];

    %%%%%%%%%%%%%%
    %% Features %%
    %%%%%%%%%%%%%%
    params.features.erp_iscompute = true;
    params.features.diffwave_iscompute = true;

    %%%%%%%%%%%%%%%%%%%%%%
    %% Resampling Ratio %%
    %%%%%%%%%%%%%%%%%%%%%%
    params.resample.is_compute = true;
    params.resample.ratio = round(params.fsamp / 64);  % re-sampling frequency is 64 Hz
    % params.resample.time = round(0.2*params.fsamp)+1:round(0.6*params.fsamp);
    params.resample.time = round(0.15*params.fsamp)+1:round(0.5*params.fsamp);
    params.resample.time = params.resample.time + params.epochOnset;

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% Power Spectral Density %%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%
    params.psd.is_compute = false;
    params.psd.type = 'pwelch';  % {'pwelch', ''pmtm', 'stockwell'}
    params.psd.time = round(0.15*params.fsamp)+1:round(0.5*params.fsamp);
    params.psd.time = params.psd.time + params.epochOnset;
    params.psd.window = hanning(length(params.psd.time));
    params.psd.nfft  = 4*params.fsamp;
    params.psd.overlap = [];
    params.psd.freq_range = [4:2:params.spectralFilter.freqs(end)];

    %%%%%%%%%%%%%%%%%%%%%%%%%
    %% Riemannien Geometry %%
    %%%%%%%%%%%%%%%%%%%%%%%%%
    params.riemann.is_compute = false;
    params.riemann.time = round(0.2*params.fsamp)+1:round(0.5*params.fsamp);
    params.riemann.time = params.riemann.time + params.epochOnset;
    params.riemann.type = 'riemann';
    params.riemann.base = [2];
    params.riemann.is_plot = true;

    %%%%%%%%%%%%%%%%
    %% Classifier %%
    %%%%%%%%%%%%%%%%
    params.classify.is_normalize = true;
    params.classify.reduction.type = 'lasso'; % {'pca', 'fisher', 'mRMR', 'lasso', 'lasso-rLDA', 'r2', 'None'}
    params.classify.type = 'linear'; % {'SVM', 'LinearSVM', 'LDA', 'diagLDA', 'diagQuadratic', 'SLR_VAR', 'L1_SLR'}
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% Asynchronous Classification %%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    params.asynchronous.waitSample = round(0.25*params.fsamp);    
end