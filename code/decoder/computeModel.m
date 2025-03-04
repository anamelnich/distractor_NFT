%function computeModel(subjectID)
subjectID='e12'
%%%%%%%%%%%%%%%%%%%%
%% Initialization %%
%%%%%%%%%%%%%%%%%%%%
clearvars -except subjectID  sessions plot_flag;
close all; clc; rng('default');
addpath(genpath('../functions'));
%%%%%%%%%%%%%%%%%%
%% Load dataset %%
%%%%%%%%%%%%%%%%%%
dataPath = [pwd '/../../data/'];
dataInfo = dir([dataPath '/' subjectID '_20*']);
data = struct();
disp(['Loading the data from ' dataInfo.name]);
[data.calibration, data.validation, data.training, data.decoding] = loadData([dataInfo.folder '/' dataInfo.name '/*']);
delete sopen.mat
%unique(data.validation.data(:,67))
%%%%%%%%%%%%%%%%%%%%%%%%
%% Set data structure %%                   
%%%%%%%%%%%%%%%%%%%%%%%%
cfg = setParams(data.validation.header);

fields = fieldnames(data);
for i = 1:numel(fields)
    fname = fields{i};
    if isempty(data.(fname))
        data = rmfield(data, fname);
        continue;
    end
    data.(fname) = preprocessDataset(data.(fname), cfg,fname);
end
fields = fieldnames(data);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Remove non-EEG channels %%                   
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
chanRemove = {'M1','M2','EOG','FP1','FP2','FPZ'};
chanIndices = find(ismember(cfg.chanLabels,chanRemove));
cfg.chanLabels(chanIndices) = [];
for i = 1:numel(fields)
    fname = fields{i};
    data.(fname).data(:, chanIndices) = [];
end
%%%%%%%%%%%%%%%%%%%%%
%% Spectral filter %%                   
%%%%%%%%%%%%%%%%%%%%%
[b,a] = butter(cfg.spectralFilter.order, cfg.spectralFilter.freqs./(cfg.fsamp/2), 'bandpass');
cfg.spectralFilter.b = b;
cfg.spectralFilter.a = a;
for i = 1:numel(fields)
    fname = fields{i};
    data.(fname).data(:,:) = filter(b,a,data.(fname).data(:,:));
end

%%%%%%%%%%%%%%
%% Epoching %%
%%%%%%%%%%%%%%
for i = 1:numel(fields)
    fname = fields{i};
    dataStruct = data.(fname);
    epochs.data = nan(length(cfg.epochSamples), length(cfg.chanlocs), length(dataStruct.index.pos));
    epochs.labels = dataStruct.index.typ;
    epochs.file_id = nan(length(dataStruct.index.typ), 1);
    epochs.file_type = cell(length(dataStruct.index.typ), 1);

    for i_trial = 1:length(dataStruct.index.pos)
        epochs.data(:, :, i_trial) = dataStruct.data(dataStruct.index.pos(i_trial) + cfg.epochSamples, :);    
        temp = find(dataStruct.index.pos(i_trial) <= dataStruct.eof, 1, 'first');
        epochs.file_id(i_trial) = temp;
    end
    data.(fname).epochs = epochs;
    data.(fname).epochs.eof = data.(fname).eof;
end 

%%%%%%%%%%%%%%%%%%%%
%% Classification %%
%%%%%%%%%%%%%%%%%%%%
trainingdata = data.training.epochs;

n_files = length(trainingdata.eof);
trainingdata.posteriors = nan(length(trainingdata.labels), 1);
for i_file = 1:n_files
    train_index = trainingdata.file_id ~= i_file; 
    test_index = trainingdata.file_id == i_file;
    decoder = computeDecoder(trainingdata.data(:, :, train_index), trainingdata.labels(train_index), cfg);
    trainingdata.posteriors(test_index) = singleClassification(decoder, trainingdata.data(:, :, test_index), trainingdata.labels(test_index),0,decoder.leftElectrodeIndices,decoder.rightElectrodeIndices);
end



%%
trainingdata.newlabels = trainingdata.labels;
trainingdata.newlabels(trainingdata.labels == 2) = 1;
disp('== Synchronous Classification == ');
% [x, y, t, auc, opt] = perfcurve(~epochs.labels,1-epochs.posteriors, 1, 'Prior', 'uniform');
[x, y, t, auc, opt] = perfcurve(trainingdata.newlabels,trainingdata.posteriors, 1, 'Prior', 'uniform');
% threshold = 0.5;
threshold = t(x == opt(1) & y == opt(2));
disp(['AUC score : ' num2str(auc, '%.2f') ' Threshold: ' num2str(threshold, '%.2f')]);
disp('Confusion Matrix: ');
cm = confusionmat(logical(trainingdata.newlabels), (trainingdata.posteriors >= threshold));
disp(cm);
tnr = cm(1,1) / sum(cm(1, :));
tpr = cm(2,2) / sum(cm(2, :));
disp(['TNR: ' num2str(tnr, '%.2f') ' TPR: ' num2str(tpr, '%.2f')]);

%%
% figure;
% histogram(epochs.posteriors(epochs.labels == 0),'FaceColor', 'green'); hold on;
% histogram(epochs.posteriors(epochs.labels == 1),'FaceColor','red'); hold on;
% histogram(epochs.posteriors(epochs.labels == 2),'FaceColor','blue');
% legend('no distractor','distractor rigt','distractorleft')
% threshold = 0.5 ;
%%
decoder = computeDecoder(trainingdata.data, trainingdata.labels, cfg);
trainingdata.posteriors = singleClassification(decoder, trainingdata.data, trainingdata.labels, 0, decoder.leftElectrodeIndices,decoder.rightElectrodeIndices);
%%
decoder.decision_threshold = threshold;
decoder.eegChannels = cfg.eegChannels; 
decoder.eogChannels = cfg.eogChannels;
%decoder.eogFilter = cfg.eogFilter;
decoder.spectralFilter = cfg.spectralFilter;
decoder.threshold = threshold;
decoder.chantoremove = chanIndices;

decoder.resample.time = decoder.resample.time - decoder.epochOnset;
decoder.psd.time = decoder.psd.time - decoder.epochOnset;
decoder.riemann.time = decoder.riemann.time - decoder.epochOnset;
decoder.asynchronous.threshold = 0.8;
decoder.asynchronous.filterLength = 1;

decoder.subjectID = subjectID;
decoder.datetime = datetime;
disp(' ');
disp('Decoder Updated at');
disp(decoder.datetime);
    
save('../cnbiLoop/decoder.mat', 'decoder');
%end
