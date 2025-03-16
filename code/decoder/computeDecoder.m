function decoder = computeDecoder(trainEpochs, trainLabels, params)
%%
rng(1)
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
%% Spatial Filter %%
%%%%%%%%%%%%%%%%%%%%
if strcmp(params.spatialFilter.type, 'CAR');
    filterMatrix = get_spatial_filter('CAR', trainEpochs,trainLabels,params);
    trainEpochs = apply_spatialFilter(trainEpochs, filterMatrix);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Baseline correction %%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
baseline_period = [-0.2, 0];
baseline_indices = find(params.epochTime >= baseline_period(1) & params.epochTime <= baseline_period(2));
baseline = mean(trainEpochs(baseline_indices, :, :), 1); % [1 x channels x trials]
trainEpochs = trainEpochs - baseline;

%%%%%%%%%%%%%%%%%%%
%% ROI selection %%
%%%%%%%%%%%%%%%%%%%

LeftElectrodes = {'P1', 'P3', 'P5', 'P7', 'PO3','PO5','PO7'};
RightElectrodes = {'P2', 'P4', 'P6', 'P8', 'PO4', 'PO6','PO8'};
% LeftElectrodes = {'PO7'};
% RightElectrodes = {'PO8'};
leftElectrodeIndices = find(ismember(params.chanLabels,LeftElectrodes));
rightElectrodeIndices = find(ismember(params.chanLabels,RightElectrodes));

n_samples = size(trainEpochs, 1);
n_trials = size(trainEpochs, 3);
n_electrodes = length(LeftElectrodes); 
selectedEpochs = nan(n_samples, n_electrodes, n_trials);

idx0 = find(trainLabels == 0);
N0 = length(idx0);
perm0 = randperm(N0);
Nleft0 = floor(N0/2);             % Number assigned to left
Nright0 = N0 - Nleft0;            % Remainder assigned to right
idxLeft0 = idx0( perm0(1:Nleft0) );
idxRight0 = idx0( perm0(Nleft0+1:end) );

for i_trial = 1:n_trials
    label = trainLabels(i_trial);

    % Select electrode indices based on the label
    if label == 1 %|| label == 0 % Distractor on right
        electrodeIndices = leftElectrodeIndices;
        % selectedEpochs(:, :, i_trial) = trainEpochs(:, leftElectrodeIndices, i_trial);
    elseif label == 2 %|| label == 0 % Distractor on left or no distractor
        electrodeIndices = rightElectrodeIndices;
        % selectedEpochs(:, :, i_trial) = trainEpochs(:, rightElectrodeIndices, i_trial);
    elseif label == 0
         if ismember(i_trial, idxLeft0)
            electrodeIndices = leftElectrodeIndices;
            % countLeft0 = countLeft0 + 1;
        else
            electrodeIndices = rightElectrodeIndices;
            % countRight0 = countRight0 + 1; 
        end
        % for e = 1:n_electrodes
        %     leftChanData  = trainEpochs(:, leftElectrodeIndices(e),  i_trial);
        %     rightChanData = trainEpochs(:, rightElectrodeIndices(e), i_trial);
        % 
        %     Take element-wise average
        %     selectedEpochs(:, e, i_trial) = 0.5 * (leftChanData + rightChanData);
        % end
    else
        error('Unknown label');
    end

    % Store the data in epochs.data
    selectedEpochs(:, :, i_trial) = trainEpochs(:, electrodeIndices, i_trial);

end


classEpochs = selectedEpochs;
%relabel epochs to be 1 (Distractor) and 0 (No Distractor)
trainLabels(trainLabels == 2) = 1;


% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %% Baseline correction %%%%%%%%%%%%%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% baseline_period = [-0.2, 0];
% baseline_indices = find(params.epochTime >= baseline_period(1) & params.epochTime <= baseline_period(2));
% baseline = mean(classEpochs(baseline_indices, :, :), 1); % [1 x channels x trials]
% classEpochs = classEpochs - baseline;

%% Spatial Filter
filterMatrix = get_spatial_filter('CCA', classEpochs, trainLabels, params);

filterMatrix = filterMatrix(:, 1:params.spatialFilter.nComp);
classEpochs = apply_spatialFilter(classEpochs, filterMatrix);
%% Temporal Information
if (params.resample.is_compute)
    resamps = classEpochs(params.resample.time(1:params.resample.ratio:end), :, :); 
    resamps = reshape(resamps, [size(resamps,1)*size(resamps,2) size(resamps,3)]);
end

%% Power Spectral Density
if (params.psd.is_compute)
    psdEpochs = classEpochs(params.psd.time, :, :);
    [psds, params] = compute_psd(params.psd.type, psdEpochs, params);
    psds = reshape(psds, [size(psds,1)*size(psds,2) size(psds,3)]);
    
    if any(isnan(psds))
        logicalIdx  = not(isnan(psds(:,1)));
        psds = psds(logicalIdx,:);
    end
end

%% Riemannien Geometry
if (params.riemann.is_compute)    
    riemannEpochs = classEpochs(params.riemann.time, :, :);
    [riemanns, params] = compute_riemann(riemannEpochs, trainLabels, params);
end

%% Concatenate all computed features
classifierEpochs = [];
startSample = 1;
if (params.resample.is_compute)
    classifierEpochs = cat(1, classifierEpochs, resamps);
    params.resample.range = startSample:size(classifierEpochs,1);
    startSample = params.resample.range(end) + 1;
end
if (params.psd.is_compute)
    classifierEpochs = cat(1, classifierEpochs, psds);
    params.psd.range = startSample:size(classifierEpochs,1);
    startSample = params.psd.range(end) + 1;
end
if (params.riemann.is_compute)
    classifierEpochs = cat(1, classifierEpochs, riemanns);
    params.riemann.range = startSample:size(classifierEpochs,1);
    startSample = params.riemann.range(end) + 1;
end

if isequal(params.classify.reduction.type, 'pca')
    [coeff, ~, ~, ~, explainedVar, mu] = pca(classifierEpochs');
    keepIdx = find(cumsum(explainedVar) > 95) - 1;
    coeff = coeff(:, 1:keepIdx);
    applyPCA = @(x) bsxfun(@minus, x', mu)*coeff;
    classifierEpochs = applyPCA(classifierEpochs)';
end

%% Normailize Features
if (params.classify.is_normalize)
    maxNorm = max(classifierEpochs, [], 2); %max value for each feature
    minNorm = min(classifierEpochs, [], 2);
    
    funNormalize = @(x) (x - minNorm) ./ (maxNorm - minNorm); %min value becomes zero, then scale to range 0 to 1
    classifierEpochs = funNormalize(classifierEpochs);
end

if isequal(params.classify.reduction.type, 'lasso')
    lambdaMax = 0.1;
    Lambda = logspace(log10(0.001*lambdaMax),log10(lambdaMax),100);
    cvmodel = fitrlinear(classifierEpochs,trainLabels,'ObservationsIn','columns', 'Lambda', Lambda, 'KFold', 5, 'Learner','leastsquares','Solver','sparsa','Regularization','lasso');
    mse = kfoldLoss(cvmodel);
    [~, idx] = min(mse);
    model = fitrlinear(classifierEpochs,trainLabels,'ObservationsIn','columns', 'Lambda', Lambda(idx), 'Learner','leastsquares','Solver','sparsa','Regularization','lasso');
    keepIdx = model.Beta~=0;
    classifierEpochs = classifierEpochs(keepIdx, :);
    
elseif isequal(params.classify.reduction.type, 'r2')
    power = compute_r2(permute(classifierEpochs, [1 3 2]), trainLabels);
    [~, keepIdx] = sort(power, 'descend');
%     keepIdx = keepIdx(1:40);
    classifierEpochs = classifierEpochs(keepIdx, :);
end

%% Classification %%
model = fitcdiscr(classifierEpochs', trainLabels, 'Prior', 'uniform', 'DiscrimType', params.classify.type);
decoder.Classes = model.ClassNames;
% q = model.Coeffs(2,1).Quadratic;
w = model.Coeffs(2,1).Linear;
mu = model.Coeffs(2,1).Const;
% distance = classifierEpochs' * q * classifierEpochs + classifierEpochs' * w + mu;
distance = classifierEpochs' * w + mu;
% 
p1 = 0.015;
p2 = 1-p1;
bcoeff1=-log((1-p1)/p1)/prctile(distance,100*p1);
bcoeff2=-log((1-p2)/p2)/prctile(distance,100*p2);
b = (bcoeff1+bcoeff2)/2;

model = @(x) 1./(1+exp(-b*(x*w+mu)));
% 
% b = 1/std(distance);
% model = @(x) 1 ./ (1 + exp(-b * (x*w + mu)));

%% Keep important functions
decoder.fsamp = params.fsamp;
decoder.epochOnset = params.epochOnset;
decoder.numFeatures = size(classifierEpochs, 1);
decoder.epochRejection = params.epochRejection;
if strcmp(params.spatialFilter.type, 'CAR');
    decoder.spatialFilter = filterMatrix;
end
decoder.spatialFilter = filterMatrix;
decoder.resample = params.resample;
decoder.psd = params.psd;
decoder.riemann = params.riemann;
decoder.classify.is_normalize = params.classify.is_normalize;
if (decoder.classify.is_normalize)
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


end

function data_output = apply_spatialFilter(data_input, filter_matrix)
[n_samples, ~, n_trials] = size(data_input);

data_output = nan(n_samples, size(filter_matrix,2), n_trials);
for i_trial = 1:n_trials
    data_output(:,:,i_trial) = data_input(:,:,i_trial) * filter_matrix;
end
end