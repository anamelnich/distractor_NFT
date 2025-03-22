function data = computeGrandAvg(trainEpochs, trainLabels, fileID, params, session)
%data.(subjectID).(fname) = computeGrandAvg(data.(fname).epochs.data, data.(fname).epochs.labels,data.(fname).epochs.file_id,cfg,[subjectID ' ' fname]);

%%%%%%%%%%%%%%%%%%%%
%% Spatial Filter %%
%%%%%%%%%%%%%%%%%%%%
if strcmp(params.spatialFilter.type, 'CAR');
    filterMatrix = get_spatial_filter('CAR', trainEpochs,trainLabels,params);
    trainEpochs = apply_spatialFilter(trainEpochs, filterMatrix);
end

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
%     fileID(badEpochIndices) = [];
% 
%     disp([num2str(sum(badEpochs)) ' / ' num2str(nTrials) ' trials are removed: ' num2str(100*sum(badEpochs)/nTrials) ' %']);
% 
% end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Baseline correction %%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
baseline_period = [-0.2, 0];
baseline_indices = find(params.epochTime >= baseline_period(1) & params.epochTime <= baseline_period(2));
baseline = mean(trainEpochs(baseline_indices, :, :), 1); % [1 x channels x trials]
trainEpochs = trainEpochs - baseline;

data.preprocEpoch = trainEpochs;
%%%%%%%%%%%%%%%%%%%
%% ROI selection %%
%%%%%%%%%%%%%%%%%%%

LeftElectrodes = {'P1', 'P3', 'P5', 'P7', 'PO3','PO5','PO7'};
RightElectrodes = {'P2', 'P4', 'P6', 'P8', 'PO4', 'PO6','PO8'};
leftElectrodeIndices = find(ismember(params.chanLabels,LeftElectrodes));
rightElectrodeIndices = find(ismember(params.chanLabels,RightElectrodes));
PO7index = find(ismember(params.chanLabels,'PO7'));
PO8index = find(ismember(params.chanLabels,'PO8'));


%%%%%%%%%%%%%%%%
%% Plot ERPs %%
%%%%%%%%%%%%%%%%
%computeGrandAvg(data.(fname).epochs.data, data.(fname).epochs.labels,data.(fname).epochs.file_id,cfg,[subjectID ' ' fname]);

path = '../../Figures/';
ROIavg = calculateGrandAvg(trainEpochs, trainLabels, leftElectrodeIndices, rightElectrodeIndices);
PO78avg = calculateGrandAvg(trainEpochs, trainLabels, PO7index, PO8index);
title1 = append(session,'all electrodes');
plotERP(ROIavg, params, title1,path)
title2 = append(session,'PO7/8');
plotERP(PO78avg,params,title2,path)

allElect_dleft = mean(trainEpochs(:, :, trainLabels == 2), 3);
allElect_dright = mean(trainEpochs(:, :, trainLabels == 1), 3);
allElect_d = mean(trainEpochs(:, :, trainLabels == 1 | trainLabels == 2), 3);
allElect_nd = mean(trainEpochs(:, :, trainLabels == 0), 3);


data.preprocEpoch = trainEpochs;
data.labels = trainLabels;
data.file_id = fileID;
data.ROIerp = ROIavg;
data.PO78erp = PO78avg;
data.allElect.allElect_d = allElect_d;
data.allElect.allElect_nd = allElect_nd;
data.allElect.allElect_dleft = allElect_dleft;
data.allElect.allElect_dright = allElect_dright;

end

function data_output = apply_spatialFilter(data_input, filter_matrix)
[n_samples, ~, n_trials] = size(data_input);

data_output = nan(n_samples, size(filter_matrix,2), n_trials);
for i_trial = 1:n_trials
    data_output(:,:,i_trial) = data_input(:,:,i_trial) * filter_matrix;
end
end