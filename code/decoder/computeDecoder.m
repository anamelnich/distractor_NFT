function decoder = computeDecoder(trainEpochs, trainLabels, params)
%%

% figure;
% plot(mean(trainEpochs(:, leftElectrodeIndices, trainLabels == 1), 3),'Color','r'); hold on;
% plot(mean(trainEpochs(:, rightElectrodeIndices, trainLabels == 1), 3),'Color','g');

% %% Spatial Filter
%remove M1 M2 EOG
%index=[1,2,3,13,19,32];
%trainEpochs(:,index,:)=[];

% filterMatrix_CAR = get_spatial_filter('CAR', trainEpochs, trainLabels, params);
% trainEpochs = apply_spatialFilter(trainEpochs, filterMatrix_CAR);
%%
% figure;
% plot(mean(trainEpochs(:, leftElectrodeIndices, trainLabels == 1), 3),'Color','r'); hold on;
% plot(mean(trainEpochs(:, rightElectrodeIndices, trainLabels == 1), 3),'Color','g');
%% Select electrodes based on epoch labels

% Initialize the array for selected electrodes
n_samples = size(trainEpochs, 1);
n_trials = size(trainEpochs, 3);
n_electrodes = 6; % Number of electrodes per trial
selectedEpochs = nan(n_samples, n_electrodes, n_trials);

% Define the electrode names for left and right
LeftElectrodes = {'P1', 'P3', 'P5', 'P7', 'PO3', 'PO5','PO7'};
RightElectrodes = {'P2', 'P4', 'P6', 'P8', 'PO4', 'PO6','PO8'};
% Map electrode names to indices
% leftElectrodeIndices = find(ismember(training.header.Label,LeftElectrodes))
% rightElectrodeIndices = find(ismember(training.header.Label,RightElectrodes))

leftElectrodeIndices = [24,25,50,51,54,62]; 
rightElectrodeIndices = [27,28,52,53,56,63]; 

% leftElectrodeIndices = [22,23,47,48,51,52,59];
% rightElectrodeIndices = [25,26,49,50,53,54,60];
% leftElectrodeIndices = [22,23,47,48,51,52,59,25,26,49,50,53,54,60];
% rightElectrodeIndices = [25,26,49,50,53,54,60,22,23,47,48,51,52,59];
%leftElectrodeIndices = [19,20,44,45,48,49,56,22,23,46,47,50,51,57];
%rightElectrodeIndices = [19,20,44,45,48,49,56,22,23,46,47,50,51,57];
% leftElectrodeIndices = [19,20,44,45,48,49,56];
% rightElectrodeIndices = [22,23,46,47,50,51,57];

%%
for i_trial = 1:n_trials
    label = trainLabels(i_trial);

    % Select electrode indices based on the label
    if label == 1 || label == 0 % Distractor on right
        electrodeIndices = leftElectrodeIndices;
    elseif label == 2 %|| label == 0 % Distractor on left or no distractor
        electrodeIndices = rightElectrodeIndices;
        
    else
        error('Unknown label');
    end

    % Store the data in epochs.data
    selectedEpochs(:, :, i_trial) = trainEpochs(:, electrodeIndices, i_trial);

end
trainEpochs = selectedEpochs;
%relabel epochs to be 1 (Distractor) and 0 (No Distractor)
trainLabels(trainLabels == 2) = 1;

% %% Spatial Filter
% filterMatrix = get_spatial_filter(params.spatialFilter.type, trainEpochs, trainLabels, params);
% 
% filterMatrix = filterMatrix(:, 1:params.spatialFilter.nComp);
% classEpochs = apply_spatialFilter(trainEpochs, filterMatrix);
classEpochs = trainEpochs;


%% Compute Channel-wise Joint Probability
if (params.epochRejection.isCompute) 
    epochRej = classEpochs(params.epochRejection.time, :, :);
    [jp, ~, params.epochRejection.distribution, params.epochRejection.data2idx] = jointprob(permute(epochRej, [2 1 3]));
    [logicalIdx, lowBound, highBound] = isoutlier(jp, 'median', 2, 'thresholdFactor', 2.75);
    
    rmIdx = any(logicalIdx, 1);
    nTrials = size(trainEpochs,3);
    rejectedLabels = trainLabels(rmIdx);
    numClass0 = sum(rejectedLabels == 0);
    numClass1 = sum(rejectedLabels == 1);
    trainEpochs(:, :, rmIdx) = [];
    trainLabels(rmIdx) = [];
    
    params.epochRejection.low = lowBound;
    params.epochRejection.high = highBound;
    
    disp([num2str(sum(rmIdx)) ' / ' num2str(nTrials) ' trials are removed: ' num2str(100*sum(rmIdx)/nTrials) ' %']);
    disp(['Rejected Class 0 ' num2str(numClass0) ' Rejected Class 1 ' num2str(numClass1)]);
    
    
%     filterMatrix = get_spatial_filter(params.spatialFilter.type, trainEpochs, trainLabels, params);
%     
%     filterMatrix = filterMatrix(:, 1:params.spatialFilter.nComp);
%     classEpochs = apply_spatialFilter(trainEpochs, filterMatrix);
    classEpochs = trainEpochs;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Baseline correction %%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
baseline_period = [-0.2, 0];
baseline_indices = find(params.epochTime >= baseline_period(1) & params.epochTime <= baseline_period(2));
baseline = mean(classEpochs(baseline_indices, :, :), 1); % [1 x channels x trials]
classEpochs = classEpochs - baseline;
%%
% figure;
% plot(mean(classEpochs(:, :, trainLabels == 1), [2,3]),'Color','r'); hold on;
% plot(mean(classEpochs(:, :, trainLabels == 0), [2,3]),'Color','g'); hold on;
% plot(mean(classEpochs(:, :, trainLabels == 2), [2,3]),'Color','b'); 
% ylim([-10 10]);
% xline(256);
% xline(564);

% CCA
%filterMatrix = get_spatial_filter(params.spatialFilter.type, classEpochs, trainLabels, params); %matrix of channels x canonical components ranked in descending order
%filterMatrix = filterMatrix(:, 1:params.spatialFilter.nComp);
%classEpochs = apply_spatialFilter(trainEpochs, filterMatrix);

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

p1 = 0.025;
p2 = 1-p1;
bcoeff1=-log((1-p1)/p1)/prctile(distance,100*p1);
bcoeff2=-log((1-p2)/p2)/prctile(distance,100*p2);
b = (bcoeff1+bcoeff2)/2;

model = @(x) 1./(1+exp(-b*(x*w+mu)));

%% Keep important functions
decoder.fsamp = params.fsamp;
decoder.epochOnset = params.epochOnset;
decoder.numFeatures = size(classifierEpochs, 1);
decoder.epochRejection = params.epochRejection;
%decoder.spatialFilter = filterMatrix;
% decoder.spatialFilter_CAR = filterMatrix_CAR;
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