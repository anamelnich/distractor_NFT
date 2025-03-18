function [posterior, epoch] = singleClassification(decoder, eeg, labels, type, leftElectrodes, rightElectrodes) %singleClassification(decoder, epochs.data(:, :, test_index))
rng(1)
%epochs.posteriors(test_index) = singleClassification(decoder, epochs.data(:, :, test_index), epochs.labels(test_index),0,decoder.leftElectrodeIndices,decoder.rightElectrodeIndices);
%ex_posterior = singleClassification(decoder, stream.eeg((first_index - round(0.2*decoder.fsamp)):end, decoder.eegChannels), label_value, 1,decoder.leftElectrodeIndices,decoder.rightElectrodeIndices); %data from first index to end of buffer
%eeg 717x64 

%if ~mislocked
%    mlock
%end
% %% Spatial Filter
%remove M1 M2 EOG Fp1 Fp2 Fpz
%index=[1,2,3,13,19,32];
%eeg(:,index,:)=[];

%CAR
% n_trials = size(eeg, 3);
% sf_eeg = nan(size(eeg,1), size(decoder.spatialFilter_CAR,2), n_trials);
% for i_trial = 1:n_trials
%     sf_eeg(:,:,i_trial) = eeg(:,:,i_trial) * decoder.spatialFilter_CAR;
% end
% eeg = sf_eeg;
if type == 1 || type == 2
    eeg(:,decoder.chantoremove)=[];
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Baseline correction %%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if type == 0
    baseline_start = decoder.epochOnset - round(0.2*decoder.fsamp); %256 - 102
    baseline = mean(eeg(baseline_start:decoder.epochOnset, :, :), 1); % [1 x channels x trials]
    eeg = eeg - baseline; 
elseif type == 1
    first_index = round(0.2*decoder.fsamp); %102
    baseline_period = [1, first_index];
    baseline = mean(eeg(baseline_period, :, :), 1); % [1 x channels x trials]
    eeg = eeg - baseline;
elseif type == 2
    eeg = eeg;
end


%% Select electrodes based on epoch labels
idx0 = find(labels == 0);
perm0 = randperm(length(idx0));
Nleft0 = floor(length(idx0)/2);
idxLeft0 = idx0(perm0(1:Nleft0));
% idxRight0 will be the complement:
idxRight0 = setdiff(idx0, idxLeft0);

n_samples = size(eeg, 1);
n_electrodes = length(leftElectrodes); 
if type == 0
    n_trials = size(eeg, 3);
    selectedEpochs = nan(n_samples, n_electrodes, n_trials);
    for i_trial = 1:n_trials
        label = labels(i_trial);
    
        % Select electrode indices based on the label
        if label == 1  %|| label == 0 % Distractor on right
            electrodeIndices = leftElectrodes;
            % selectedEpochs(:, :, i_trial) = eeg(:, leftElectrodes, i_trial);
        elseif label == 2 %|| label == 0 % Distractor on left or no distractor
            electrodeIndices = rightElectrodes;
            % selectedEpochs(:, :, i_trial) = eeg(:, rightElectrodes, i_trial);
        elseif label == 0
            if ismember(i_trial, idxLeft0)
                electrodeIndices = leftElectrodes;
                % countLeft0 = countLeft0 + 1;
            else
                electrodeIndices = rightElectrodes;
                % countRight0 = countRight0 + 1; 
            end
            % for e = 1:n_electrodes
            %     leftChanData  = eeg(:, leftElectrodes(e),  i_trial);
            %     rightChanData = eeg(:, rightElectrodes(e), i_trial);
            % 
            %     % Take element-wise average
            %     selectedEpochs(:, e, i_trial) = 0.5 * (leftChanData + rightChanData);
            % end
        else
            error('Unknown label');
        end
    
        % Store the data in epochs.data
        selectedEpochs(:, :, i_trial) = eeg(:, electrodeIndices, i_trial);
    
    end
    labels(labels == 2) = 1;
elseif type == 1
    selectedEpochs = nan(n_samples, n_electrodes);
    n_trials = 1;

    % Select electrode indices based on the label
    if labels == 102
        electrodeIndices = leftElectrodes;
    elseif labels == 104 %|| label == 0 % Distractor on left or no distractor
        electrodeIndices = rightElectrodes;
    % elseif labels == 100 || labels == 110
    
    else
        error('Unknown label');
    end

    % Store the data in epochs.data
    selectedEpochs(:, :) = eeg(:, electrodeIndices); %717x7
elseif type == 2
    selectedEpochs(:,:) = eeg(:,rightElectrodes);
    n_trials = 1;
end
% disp(['Label=0: ' num2str(length(idxLeft0)) ' trials used LEFT electrodes, ' ...
%       num2str(length(idxRight0)) ' trials used RIGHT electrodes.']);
eeg = selectedEpochs;


% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %% Baseline correction %%%%%%%%%%%%%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% if type == 0
%     baseline_start = decoder.epochOnset - round(0.2*decoder.fsamp); %256 - 102
%     baseline = mean(eeg(baseline_start:decoder.epochOnset, :, :), 1); % [1 x channels x trials]
%     eeg = eeg - baseline; 
% elseif type == 1
%     first_index = round(0.2*decoder.fsamp); %102
%     baseline_period = [1, first_index];
%     baseline = mean(eeg(baseline_period, :, :), 1); % [1 x channels x trials]
%     eeg = eeg - baseline;
%    % epochStart = first_index - round(0.5*decoder.fsamp);
%     %eeg = eeg(epochStart:end,:,:);
%     %decoder.resample.time = round(0.1*params.fsamp)+1:round(0.6*params.fsamp); %check if this makes sense to do
% %     sf_eeg = eeg(first_index:end,:,:);
% elseif type == 2
%     eeg = eeg;
% %     sf_eeg = eeg;
% end
% % sf_eeg = eeg;

%% Spatial Filter
n_trials = size(eeg, 3);
sf_eeg = nan(size(eeg,1), size(decoder.spatialFilter,2), n_trials);
for i_trial = 1:n_trials
    sf_eeg(:,:,i_trial) = eeg(:,:,i_trial) * decoder.spatialFilter;
end

%% Temporal Information
if (decoder.resample.is_compute)
    if type == 0 || type == 2
        resamp = sf_eeg(decoder.resample.time(1:decoder.resample.ratio:end), :, :); %334 to 512 or ~0.65-0.5 = 0.15 to 1-0.5=0.5
        resamp = reshape(resamp, [size(resamp,1)*size(resamp,2) n_trials]);
    elseif type == 1
        time = first_index + (round(0.15*decoder.fsamp)+1:round(0.5*decoder.fsamp)); %78 through 358
        resamp = sf_eeg(time(1:decoder.resample.ratio:end), :, :); %this doesn't seem right for online cause resample.time(1) = 308
        resamp = reshape(resamp, [size(resamp,1)*size(resamp,2) n_trials]);
    end
end

%% Power Spectral Density
if (decoder.psd.is_compute)
    psd_epoch = sf_eeg(decoder.psd.time, :,:);
    [psd, decoder] = compute_psd(decoder.psd.type, psd_epoch, decoder);
    psd = reshape(psd, [size(psd,1)*size(psd,2) n_trials]);
end

%% Riemannien Geometry
if (decoder.riemann.is_compute)
    riemann_epoch = sf_eeg(decoder.riemann.time, :, :);
    cov_matrix = covariances(permute(cat(2, repmat(decoder.riemann.avrgSignals, [1 1 n_trials]), riemann_epoch), [2 1 3]), 'shcovft');
    riemann = Tangent_space(cov_matrix, decoder.riemann.avrgCov);
end

%% Concatenate all computed features
epoch = nan(decoder.numFeatures, n_trials);
if (decoder.resample.is_compute)
    epoch(decoder.resample.range,:) = resamp;
end
if (decoder.psd.is_compute)
    epoch(decoder.psd.range,:) = psd;
end
if (decoder.riemann.is_compute)
    epoch(decoder.riemann.range,:) = riemann;
end

if isequal(decoder.classify.reduction.type, 'pca')
    epoch = decoder.classify.applyPCA(epoch)';
end

%% Normailize Features
if (decoder.classify.is_normalize)
    epoch = decoder.classify.funNormalize(epoch);
end

if ismember(decoder.classify.reduction.type, {'lasso', 'r2'})
    epoch = epoch(decoder.classify.keepIdx, :);
end

%% Classification
posterior = decoder.classify.model(epoch');
