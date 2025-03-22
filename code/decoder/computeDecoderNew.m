function [decoder, classifierEpochs] = computeDecoderNew(trainEpochs, trainLabels, params)
% computeDecoderNew(trainingData.data, trainingData.labels, cfg);
%%
% rng(1)

%%%%%%%%%%%%%%%%%%%%%%%%
%% Artifact Rejection %%
%%%%%%%%%%%%%%%%%%%%%%%%
% if (params.epochRejection.isCompute) 
%     badEpochs = artRej(trainEpochs, params);
% 
%     badEpochIndices = find(badEpochs);
%     nTrials = size(trainEpochs,3);
%     trainEpochs(:, :, badEpochIndices) = [];
%     trainLabels(badEpochIndices) = [];
% 
% 
%     disp([num2str(sum(badEpochs)) ' / ' num2str(nTrials) ' trials are removed: ' num2str(100*sum(badEpochs)/nTrials) ' %']);
% 
% end

%%%%%%%%%%%%%%%%%%%%
%% Spatial Filter (CAR) %%
%%%%%%%%%%%%%%%%%%%%
if strcmp(params.spatialFilter.type, 'CAR')
    filterMatrix = get_spatial_filter('CAR', trainEpochs, trainLabels, params);
    trainEpochs = apply_spatialFilter(trainEpochs, filterMatrix);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Baseline Correction %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
baseline_period = [-0.2, 0];
baseline_indices = find(params.epochTime >= baseline_period(1) & params.epochTime <= baseline_period(2));
baseline = mean(trainEpochs(baseline_indices, :, :), 1);  % [1 x channels x trials]
trainEpochs = trainEpochs - baseline;

%%%%%%%%%%%%%%%%%%%
%% ROI Selection %%
%%%%%%%%%%%%%%%%%%%
LeftElectrodes = {'P1', 'P3', 'P5', 'P7', 'PO3', 'PO5', 'PO7'};
RightElectrodes = {'P2', 'P4', 'P6', 'P8', 'PO4', 'PO6', 'PO8'};
% LeftElectrodes = {'P1', 'P3', 'P5', 'P7', 'PO3', 'PO5', 'PO7','P2', 'P4', 'P6', 'P8', 'PO4', 'PO6', 'PO8'};
% RightElectrodes = {'P1', 'P3', 'P5', 'P7', 'PO3', 'PO5', 'PO7','P2', 'P4', 'P6', 'P8', 'PO4', 'PO6', 'PO8'};
leftElectrodeIndices = find(ismember(params.chanLabels, LeftElectrodes));
rightElectrodeIndices = find(ismember(params.chanLabels, RightElectrodes));

n_samples = size(trainEpochs, 1);
n_trials = size(trainEpochs, 3);
n_electrodes = length(LeftElectrodes);

% Pre-calculate indices for label==0 random assignment (for ERP & diff)
idx0 = find(trainLabels == 0);
perm0 = randperm(length(idx0));
Nleft0 = floor(length(idx0)/2);
idxLeft0 = idx0(perm0(1:Nleft0));
% idxRight0 will be the complement:
idxRight0 = setdiff(idx0, idxLeft0);


% Allocate matrices for ERP and difference wave
selectedEpochsERP = nan(n_samples, n_electrodes, n_trials);
selectedEpochsDiff = nan(n_samples, n_electrodes, n_trials);

for i_trial = 1:n_trials
    label = trainLabels(i_trial);
    
    % --- ERP (regular) selection ---
    if label == 1
        electrodeIndicesERP = leftElectrodeIndices;
    elseif label == 2
        electrodeIndicesERP = rightElectrodeIndices;
    elseif label == 0
        if rand < 0.5
%         if ismember(i_trial, idxLeft0)
            electrodeIndicesERP = leftElectrodeIndices;
        else
            electrodeIndicesERP = rightElectrodeIndices;
        end
    else
        error('Unknown label');
    end
    selectedEpochsERP(:, :, i_trial) = trainEpochs(:, electrodeIndicesERP, i_trial);
    
    % --- Difference wave computation ---
    % For each trial, subtract corresponding electrodes.
    if label == 1
        diff_trial = trainEpochs(:, leftElectrodeIndices, i_trial) - trainEpochs(:, rightElectrodeIndices, i_trial);
    elseif label == 2
        diff_trial = trainEpochs(:, rightElectrodeIndices, i_trial) - trainEpochs(:, leftElectrodeIndices, i_trial);
    elseif label == 0
        if rand < 0.5
%         if ismember(i_trial, idxLeft0)
            diff_trial = trainEpochs(:, leftElectrodeIndices, i_trial) - trainEpochs(:, rightElectrodeIndices, i_trial);
        else
            diff_trial = trainEpochs(:, rightElectrodeIndices, i_trial) - trainEpochs(:, leftElectrodeIndices, i_trial);
        end
    end
    selectedEpochsDiff(:, :, i_trial) = diff_trial;
end
trainLabels(trainLabels == 2) = 1;
% At this point, selectedEpochsERP and selectedEpochsDiff are both
% [n_samples x n_electrodes x n_trials].

% Process ERP features if requested.
if params.features.erp_iscompute
    [ERP_classEpochs, ERPfilterMatrix] = processFeatures(selectedEpochsERP, trainLabels, params); % 46x240
else
    ERP_classEpochs = [];
    ERPfilterMatrix = 'na';
end

% Process difference-wave features if requested.
if params.features.diffwave_iscompute
    [diff_classEpochs, diffFilterMatrix] = processFeatures(selectedEpochsDiff, trainLabels, params);%46x240
else
    diff_classEpochs = [];
    diffFilterMatrix = 'na';
end

% Concatenate along feature dimension.
if params.features.erp_iscompute && params.features.diffwave_iscompute
    classifierEpochs = [ERP_classEpochs; diff_classEpochs]; %92x240
elseif params.features.erp_iscompute
    classifierEpochs = ERP_classEpochs;
elseif params.features.diffwave_iscompute
    classifierEpochs = diff_classEpochs;
else
    error('No features selected. Set either params.features.erp or params.features.diffwave to true.');
end


if isequal(params.classify.reduction.type, 'pca')
    [coeff, ~, ~, ~, explainedVar, mu] = pca(classifierEpochs');
    keepIdx = find(cumsum(explainedVar) > 95) - 1;
    coeff = coeff(:, 1:keepIdx);
    applyPCA = @(x) bsxfun(@minus, x', mu)*coeff;
    classifierEpochs = applyPCA(classifierEpochs)';
end

if (params.classify.is_normalize)
    maxNorm = max(classifierEpochs, [], 2); % max value for each feature
    minNorm = min(classifierEpochs, [], 2);
    funNormalize = @(x) (x - minNorm) ./ (maxNorm - minNorm);
    classifierEpochs = funNormalize(classifierEpochs);
end



if isequal(params.classify.reduction.type, 'lasso')
    lambdaMax = 0.1;
    Lambda = logspace(log10(0.001*lambdaMax), log10(lambdaMax), 100);
    cvmodel = fitrlinear(classifierEpochs, trainLabels, 'ObservationsIn', 'columns', 'Lambda', Lambda, 'KFold', 5, 'Learner', 'leastsquares', 'Solver', 'sparsa', 'Regularization', 'lasso');
    mse = kfoldLoss(cvmodel);
    [~, idx] = min(mse);
    selectedLambda = Lambda(idx);
    model = fitrlinear(classifierEpochs, trainLabels, 'ObservationsIn', 'columns', 'Lambda', Lambda(idx), 'Learner', 'leastsquares', 'Solver', 'sparsa', 'Regularization', 'lasso');
    keepIdx = model.Beta ~= 0;
    classifierEpochs = classifierEpochs(keepIdx, :);
    disp(['Number of features selected: ', num2str(sum(keepIdx))]);
elseif isequal(params.classify.reduction.type, 'r2')
    power = compute_r2(permute(classifierEpochs, [1 3 2]), trainLabels); %after permute, 92x1x240, trainLabels 240x1
    % [~, keepIdx] = sort(power, 'descend');
    % keepIdx = keepIdx(1:30);

    maxPower = max(power);
    disp(maxPower)
    threshold = 0.25 * maxPower;
    disp(threshold)
    keepIdx = find(power >= threshold);
    disp(length(keepIdx))

    classifierEpochs = classifierEpochs(keepIdx, :);
    
    figure;
    bar(1:length(power), power);
    hold on;
    bar(keepIdx, power(keepIdx));
    xlabel('Feature Index');
    ylabel('r^2');
    ylim([0 0.1]);
    title('Feature-wise r^2 (with selected features highlighted)');
    hold off;
    figpath = '../../Figures/';
    baseName = 'Featurewise_r2.png';
    fullName = fullfile(figpath, baseName);
    counter = 1;
    while exist(fullName, 'file')
        fullName = fullfile(figpath, sprintf('Featurewise_r2_%d.png', counter));
        counter = counter + 1;
    end
    saveas(gcf, fullName);
end

%% Classification %%
% disp(['Size of classifierEpochs: ', mat2str(size(classifierEpochs))]);
% disp(['Size of trainLabels: ', mat2str(size(trainLabels))]);

model = fitcdiscr(classifierEpochs', trainLabels, 'Prior', 'uniform', 'DiscrimType', params.classify.type);
decoder.Classes = model.ClassNames;
w = model.Coeffs(2,1).Linear;
mu_coef = model.Coeffs(2,1).Const;
distance = classifierEpochs' * w + mu_coef;
p1 = 0.025;
p2 = 1-p1;
bcoeff1 = -log((1-p1)/p1)/prctile(distance, 100*p1);
bcoeff2 = -log((1-p2)/p2)/prctile(distance, 100*p2);
b = (bcoeff1+bcoeff2)/2;
model = @(x) 1./(1+exp(-b*(x*w+mu_coef)));

%% Keep important parameters in decoder
decoder.fsamp = params.fsamp;
decoder.epochOnset = params.epochOnset;
decoder.numFeatures = size(classifierEpochs, 1);
decoder.epochRejection = params.epochRejection;
if strcmp(params.spatialFilter.type, 'CAR')
    decoder.spatialFilter = filterMatrix;
end
decoder.spatialFilter.erp = ERPfilterMatrix;
decoder.spatialFilter.diff = diffFilterMatrix;
decoder.resample = params.resample;
decoder.psd = params.psd;
decoder.riemann = params.riemann;
decoder.classify.is_normalize = params.classify.is_normalize;
if decoder.classify.is_normalize
    decoder.classify.funNormalize = funNormalize;
end
decoder.classify.reduction.type = params.classify.reduction.type;
if ismember(decoder.classify.reduction.type, {'lasso', 'r2'})
    decoder.classify.keepIdx = keepIdx;
elseif isequal(decoder.classify.reduction.type, 'pca')
    decoder.classify.applyPCA = applyPCA;
end
decoder.classify.type = params.classify.type;
decoder.classify.model = model;
decoder.leftElectrodeIndices = leftElectrodeIndices;
decoder.rightElectrodeIndices = rightElectrodeIndices;
decoder.features.erp_iscompute = params.features.erp_iscompute;
decoder.features.diffwave_iscompute = params.features.diffwave_iscompute;
decoder.baseline_indices = baseline_indices;
if isequal(params.classify.reduction.type, 'lasso')
    decoder.lassoLambda = selectedLambda;
end 

end

%% Helper Function: processFeatures
function [features, filterMatrix] = processFeatures(epochData, trainLabels, params)
% processFeatures applies the same preprocessing steps (spatial filtering, 
% resampling, PSD/riemann if applicable, dimensionality reduction, normalization)
% to the input epoch data.
%
% Input:
%   epochData  - a 3D matrix [time x channels x trials]
%   trainLabels- vector of labels
%   params     - configuration structure containing all settings.
%
% Output:
%   features   - a feature matrix, e.g., [numFeatures x nTrials]

% Apply CCA-based spatial filter (or any specified filter).
filterMatrix = get_spatial_filter('CCA', epochData, trainLabels, params);
filterMatrix = filterMatrix(:, 1:params.spatialFilter.nComp);
epochData = apply_spatialFilter(epochData, filterMatrix);

% Temporal Information: Resampling
if (params.resample.is_compute)
    resamps = epochData(params.resample.time(1:params.resample.ratio:end), :, :); 
    resamps = reshape(resamps, [size(resamps,1)*size(resamps,2) size(resamps,3)]);
else
    resamps = [];
end

% Power Spectral Density computation (if enabled)
if (params.psd.is_compute)
    psdEpochs = epochData(params.psd.time, :, :);
    [psds, params] = compute_psd(params.psd.type, psdEpochs, params);
    psds = reshape(psds, [size(psds,1)*size(psds,2) size(psds,3)]);
    if any(isnan(psds))
        logicalIdx = ~isnan(psds(:,1));
        psds = psds(logicalIdx, :);
    end
else
    psds = [];
end

% Riemannian Geometry features (if enabled)
if (params.riemann.is_compute)
    riemannEpochs = epochData(params.riemann.time, :, :);
    [riemanns, params] = compute_riemann(riemannEpochs, trainLabels, params);
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
