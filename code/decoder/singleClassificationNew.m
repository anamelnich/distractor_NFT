function [posterior, epoch, NDside] = singleClassificationNew(decoder, eeg, labels, type, leftElectrodes, rightElectrodes)
% singleClassification applies the same preprocessing steps as computeDecoder
% to test data and returns the posterior probability for the trial(s),
% along with the final feature vector(s) (epoch).
%
% It now also computes a difference wave:
%   - If label==1: difference = left electrodes - right electrodes.
%   - If label==2: difference = right electrodes - left electrodes.
%   - If label==0: a random (50/50) assignment is used (consistent with training).
%
% The same preprocessing (baseline correction, spatial filtering, resampling,
% PSD/riemann, PCA, normalization, etc.) is applied to both ERP and difference
% wave data. Final features are concatenated according to the flags:
%   decoder.features.erp and decoder.features.diffwave.
%
% Inputs:
%   decoder       - structure containing all transformation parameters.
%   eeg           - test EEG data. For type==0: [samples x channels x trials];
%                   for type==1 or 2: [samples x channels] (a single trial).
%   labels        - labels for the epoch(s) (vector for type==0, scalar for type==1/2).
%   type          - integer (0: offline multi-trial; 1: online; 2: special case)
%   leftElectrodes, rightElectrodes - indices or cell arrays specifying ROI electrodes.
%
% Outputs:
%   posterior     - posterior probability computed by the classifier.
%   epoch         - final feature vector(s) used by the classifier.

rng(1)  % Ensure reproducible random choices

%% (Optional) Remove unwanted channels for type==1 or type==2
if type == 1 || type == 2
    eeg(:, decoder.chantoremove) = [];
end

%% Baseline Correction
if type == 0
    baseline_start = decoder.epochOnset - round(0.2 * decoder.fsamp);
    baseline = mean(eeg(baseline_start:decoder.epochOnset, :, :), 1); % [1 x channels x trials]
    eeg = eeg - baseline;
    first_index = 0; %placeholder
elseif type == 1
    first_index = round(0.2 * decoder.fsamp);
    baseline_period = [1, first_index];
    baseline = mean(eeg(baseline_period, :, :), 1);
    eeg = eeg - baseline;
elseif type == 2
    % No baseline correction applied.
    first_index = 0; %placeholder
end

%% ROI Selection & Difference Wave Computation
% (We compute both ERP and diff wave epochs.)
NDside=100
if type == 0
    % Offline multi-trial: eeg is 3D.
    n_trials = size(eeg, 3);
    n_samples = size(eeg, 1);
    n_electrodes = length(leftElectrodes);
    selectedEpochsERP = nan(n_samples, n_electrodes, n_trials);
    selectedEpochsDiff = nan(n_samples, n_electrodes, n_trials);
    
    % Precompute random assignment for label==0:
    idx0 = find(labels == 0);
    perm0 = randperm(length(idx0));
    Nleft0 = floor(length(idx0)/2);
    idxLeft0 = idx0(perm0(1:Nleft0));
    
    for i_trial = 1:n_trials
        lab = labels(i_trial);
        % ERP: Select electrodes based on label.
        if lab == 1
            electrodeIndicesERP = leftElectrodes;
            diff_trial = eeg(:, leftElectrodes, i_trial) - eeg(:, rightElectrodes, i_trial);
        elseif lab == 2
            electrodeIndicesERP = rightElectrodes;
            diff_trial = eeg(:, rightElectrodes, i_trial) - eeg(:, leftElectrodes, i_trial);
        elseif lab == 0
            if ismember(i_trial, idxLeft0)
                electrodeIndicesERP = leftElectrodes;
                diff_trial = eeg(:, leftElectrodes, i_trial) - eeg(:, rightElectrodes, i_trial);
            else
                electrodeIndicesERP = rightElectrodes;
                diff_trial = eeg(:, rightElectrodes, i_trial) - eeg(:, leftElectrodes, i_trial);
            end
        else
            error('Unknown label');
        end
        selectedEpochsERP(:, :, i_trial) = eeg(:, electrodeIndicesERP, i_trial);
        selectedEpochsDiff(:, :, i_trial) = diff_trial;
    end
elseif type == 1
    % Online: eeg is 2D.
    n_samples = size(eeg, 1);
    n_electrodes = length(leftElectrodes);
    % ERP selection for the single trial.
    if labels == 102  % for example, meaning distractor on right -> use left electrodes
        electrodeIndicesERP = leftElectrodes;
        diff_trial = eeg(:, leftElectrodes) - eeg(:, rightElectrodes);
    elseif labels == 104
        electrodeIndicesERP = rightElectrodes;
        diff_trial = eeg(:, rightElectrodes) - eeg(:, leftElectrodes);
    elseif labels == 100 || labels == 110
        if rand < 0.5
            electrodeIndicesERP = leftElectrodes;
            diff_trial = eeg(:, leftElectrodes) - eeg(:, rightElectrodes);
            NDside = 0;
        else
            electrodeIndicesERP = rightElectrodes;
            diff_trial = eeg(:, rightElectrodes) - eeg(:, leftElectrodes);
            NDside = 1;
        end
    else
        error('Unknown label');
    end
    selectedEpochsERP = eeg(:, electrodeIndicesERP);
    selectedEpochsDiff = diff_trial;
    n_trials = 1;
elseif type == 2
    % For type 2, assume we only use ERP from right electrodes.
    n_samples = size(eeg, 1);
    selectedEpochsERP = eeg(:, rightElectrodes);
    selectedEpochsDiff = [];  % Not computed in this mode.
    n_trials = 1;
end

%% Feature Processing
% Use the same helper function (processFeatures) as in computeDecoder.
% This function applies spatial filtering (e.g., using CCA), resampling, PSD/riemann, etc.
if decoder.features.erp_iscompute
    ERP_features = processFeatures(selectedEpochsERP, decoder, decoder.spatialFilter.erp,first_index);
else
    ERP_features = [];
end
if decoder.features.diffwave_iscompute
    diff_features = processFeatures(selectedEpochsDiff, decoder, decoder.spatialFilter.diff,first_index);
else
    diff_features = [];
end

if ~isempty(ERP_features) && ~isempty(diff_features)
    epoch = [ERP_features; diff_features];
elseif ~isempty(ERP_features)
    epoch = ERP_features;
elseif ~isempty(diff_features)
    epoch = diff_features;
else
    error('No features selected in decoder.features');
end

%% Further Processing: PCA, Normalization, etc.
if isequal(decoder.classify.reduction.type, 'pca')
    epoch = decoder.classify.applyPCA(epoch)';
end
if decoder.classify.is_normalize
    epoch = decoder.classify.funNormalize(epoch);
end
if ismember(decoder.classify.reduction.type, {'lasso','r2'})
    epoch = epoch(decoder.classify.keepIdx, :);
end

%% Classification
posterior = decoder.classify.model(epoch');

end

%% Helper Function: processFeatures
function features = processFeatures(eeg, decoder, filterMatrix,type, first_index)
% Spatial Filter
n_trials = size(eeg, 3);
% Allocate a new variable for the spatially filtered EEG
sf_eeg = nan(size(eeg,1), size(filterMatrix,2), n_trials);
for i_trial = 1:n_trials
    % Multiply the trial's data [n_samples x n_channels] by filterMatrix [n_channels x nComp]
    sf_eeg(:,:,i_trial) = eeg(:,:,i_trial) * filterMatrix;
end

%% Temporal Information
if (decoder.resample.is_compute)
    if type == 0 || type == 2
        resamp = sf_eeg(decoder.resample.time(1:decoder.resample.ratio:end), :, :); %334 to 512 or ~0.65-0.5 = 0.15 to 1-0.5=0.5
        resamps = reshape(resamp, [size(resamp,1)*size(resamp,2) n_trials]);
    elseif type == 1
        time = first_index + (round(0.15*decoder.fsamp)+1:round(0.5*decoder.fsamp)); %78 through 358
        resamp = sf_eeg(time(1:decoder.resample.ratio:end), :, :); 
        resamps = reshape(resamp, [size(resamp,1)*size(resamp,2) n_trials]);
    end
else
    resamps = [];
end

% Power Spectral Density computation (if enabled)
if (decoder.psd.is_compute)
    psd_epoch = sf_eeg(decoder.psd.time, :,:);
    [psd, decoder] = compute_psd(decoder.psd.type, psd_epoch, decoder);
    psds = reshape(psd, [size(psd,1)*size(psd,2) n_trials]);
else
    psds = [];
end

% Riemannian Geometry features (if enabled)
if (decoder.riemann.is_compute)
    riemann_epoch = sf_eeg(decoder.riemann.time, :, :);
    cov_matrix = covariances(permute(cat(2, repmat(decoder.riemann.avrgSignals, [1 1 n_trials]), riemann_epoch), [2 1 3]), 'shcovft');
    riemanns = Tangent_space(cov_matrix, decoder.riemann.avrgCov);
else
    riemanns = [];
end

% Concatenate all computed features
features = [];
startSample = 1;
if ~isempty(resamps)
    features = cat(1, features, resamps);
    params.resample.range = startSample:size(features,1);
    startSample = params.resample.range(end) + 1;
end
if ~isempty(psds)
    features = cat(1, features, psds);
    params.psd.range = startSample:size(features,1);
    startSample = params.psd.range(end) + 1;
end
if ~isempty(riemanns)
    features = cat(1, features, riemanns);
    params.riemann.range = startSample:size(features,1);
    startSample = params.riemann.range(end) + 1;
end

end

%% apply_spatialFilter
function data_output = apply_spatialFilter(data_input, filter_matrix)
[n_samples, ~, n_trials] = size(data_input);
data_output = nan(n_samples, size(filter_matrix,2), n_trials);
for i_trial = 1:n_trials
    data_output(:,:,i_trial) = data_input(:,:,i_trial) * filter_matrix;
end
end