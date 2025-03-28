%function computeModel(subjectID)
subjectID='e13'
%%%%%%%%%%%%%%%%%%%%
%% Initialization %%
%%%%%%%%%%%%%%%%%%%%
clearvars -except subjectID  sessions plot_flag;
close all; clc; rng('default');
addpath(genpath('../functions'));
figpath = '../../Figures/';
%%%%%%%%%%%%%%%%%%
%% Load dataset %%
%%%%%%%%%%%%%%%%%%
dataPath = [pwd '/../../data/'];
data = loadData(dataPath, subjectID);
delete sopen.mat

%%%%%%%%%%%%%%%%%%%%%%%%
%% Set data structure %%                   
%%%%%%%%%%%%%%%%%%%%%%%%
cfg = setParams(data.eogcalibration1.header);

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

%session 1 model
epochsForTrain = {data.training1.epochs};
% epochsForTrain = {data.training1.epochs, data.decoding1.epochs};
%session 2 model
% epochsForTrain = {data.training1.epochs,data.training2.epochs,  ...
% data.decoding1.epochs};   

%session 3 model
% epochsForTrain = {data.training1.epochs, data.training2.epochs, ...
%  data.decoding1.epochs,data.decoding2.epochs};
%%
trainingData = combineEpochs(epochsForTrain);
cfg.features.diffwave_iscompute = true;
cfg.spatialFilter.nComp = 2;
cfg.classify.reduction.type = 'r2';
n_files = length(trainingData.eof);
trainingData.posteriors = nan(length(trainingData.labels), 1);
trainingData.posteriorstrain = nan(length(trainingData.labels), 1);
for i_file = 1:n_files
    train_index = trainingData.file_id ~= i_file; 
    test_index = trainingData.file_id == i_file;
    decoder = computeDecoderNew(trainingData.data(:, :, train_index), trainingData.labels(train_index), cfg);
    trainingData.posteriors(test_index) = singleClassificationNew(decoder, trainingData.data(:, :, test_index), trainingData.labels(test_index),0,decoder.leftElectrodeIndices,decoder.rightElectrodeIndices);
end



%%
trainingData.newlabels = trainingData.labels;
trainingData.newlabels(trainingData.labels == 2) = 1;
disp('== Synchronous Classification == ');
% [x, y, t, auc, opt] = perfcurve(~trainingData.newlabels,1-trainingData.posteriors, 1, 'Prior', 'uniform');
[x, y, t, auc, opt] = perfcurve(trainingData.newlabels,trainingData.posteriors, 1, 'Prior', 'uniform');
threshold = 0.4961;
% threshold = t(x == opt(1) & y == opt(2));
disp(['AUC score : ' num2str(auc, '%.2f') ' Threshold: ' num2str(threshold, '%.2f')]);
disp('Confusion Matrix: ');
cmCV = confusionmat(logical(trainingData.newlabels), (trainingData.posteriors >= threshold));
disp('Confusion Matrix (with labels):')
disp('--------------------------------')
disp('            Pred=0    Pred=1')
fprintf('True=0:       %3d       %3d\n', cmCV(1,1), cmCV(1,2));
fprintf('True=1:       %3d       %3d\n', cmCV(2,1), cmCV(2,2));
tnr = cmCV(1,1) / sum(cmCV(1, :));
tpr = cmCV(2,2) / sum(cmCV(2, :));
accuracy = (cmCV(1,1) + cmCV(2,2)) / sum(cmCV(:));
disp(['TNR: ' num2str(tnr, '%.2f') ' TPR: ' num2str(tpr, '%.2f') ' Accuracy: ' num2str(accuracy, '%.2f')]);


%%
figure;

% Subplot 1: Label = 1
subplot(4,1,1);
histogram(trainingData.posteriors(trainingData.labels == 1), ...
    'FaceAlpha', 0.5, ...
    'FaceColor', 'r', ...
    'EdgeColor', 'none', ...
    'NumBins', 10);
title('Label = 1');
xlabel('Posterior Probability');
ylabel('Count');
xlim([0 1]);
ylim([0 25]);

% Subplot 2: Label = 2
subplot(4,1,2);
histogram(trainingData.posteriors(trainingData.labels == 2), ...
    'FaceAlpha', 0.5, ...
    'FaceColor', 'g', ...
    'EdgeColor', 'none', ...
    'NumBins', 10);
title('Label = 2');
xlabel('Posterior Probability');
ylabel('Count');
xlim([0 1]);
ylim([0 25]);


% Subplot 3: Label = 0
subplot(4,1,3);
histogram(trainingData.posteriors(trainingData.labels == 0), ...
    'FaceAlpha', 0.5, ...
    'FaceColor', 'b', ...
    'EdgeColor', 'none', ...
    'NumBins', 10);
title('Label = 0');
xlabel('Posterior Probability');
ylabel('Count');
xlim([0 1]);
% ylim([0 150]);
% ylim([0 90]);
ylim([0 50]);


% Subplot 3: Label = 0
subplot(4,1,4);
histogram(trainingData.posteriors(trainingData.newlabels == 1), ...
    'FaceAlpha', 0.5, ...
    'FaceColor', 'r', ...
    'EdgeColor', 'none', ...
    'NumBins', 10);
title('Label = 1');
xlabel('Posterior Probability');
ylabel('Count');
xlim([0 1]);
% ylim([0 150]);
% ylim([0 90]);
ylim([0 50]);

%% Visualize CCA
Electrodes = {'P1/2', 'P3/4', 'P5/6', 'P7/8', 'PO3/4', 'PO5/6', 'PO7/8'};
figure;
subplot(1,2,1);
imagesc(decoder.spatialFilter.erp); % Displays a heatmap
colorbar;
title('ERP Spatial Filter Weights');
xlabel('Component');
set(gca, 'XTick', 1:2, 'XTickLabel', {'Comp 1', 'Comp 2'});
ylabel('Left Electrodes');
set(gca, 'YTick', 1:7, 'YTickLabel', Electrodes);

subplot(1,2,2);
imagesc(decoder.spatialFilter.diff);
colorbar;
title('Diff Spatial Filter Weights');
xlabel('Component');
set(gca, 'XTick', 1:2, 'XTickLabel', {'Comp 1', 'Comp 2'});
ylabel('Right Electrodes');
set(gca, 'YTick', 1:7, 'YTickLabel', Electrodes);


%%
[decoder,modeloutput] = computeDecoderNew(trainingData.data, trainingData.labels, cfg);
[trainingData.posteriors, classoutput] = singleClassificationNew(decoder, trainingData.data, trainingData.labels, 0, decoder.leftElectrodeIndices,decoder.rightElectrodeIndices);
if isequal(modeloutput, classoutput)
    disp('Preprocessing is the same.');
else
    disp('Preprocessing is NOT the same.');
end
cm = confusionmat(logical(trainingData.newlabels), (trainingData.posteriors >= threshold));
disp('Confusion Matrix (with labels):')
disp('--------------------------------')
disp('            Pred=0    Pred=1')
fprintf('True=0:       %3d       %3d\n', cm(1,1), cm(1,2));
fprintf('True=1:       %3d       %3d\n', cm(2,1), cm(2,2));
tnr = cm(1,1) / sum(cm(1, :));
tpr = cm(2,2) / sum(cm(2, :));
accuracy = (cm(1,1) + cm(2,2)) / sum(cm(:));
disp(['TNR: ' num2str(tnr, '%.2f') ' TPR: ' num2str(tpr, '%.2f') ' Accuracy: ' num2str(accuracy, '%.2f')]);

%% Offline posteriors 
%% Set Parameters
session = 3;         % Choose session: 1, 2, or 3

switch session
    case 1
        epochsForTrain = {data.training1.epochs};
        testingData = data.decoding1.epochs;
    case 2
        epochsForTrain = {data.training1.epochs, data.training2.epochs, data.decoding1.epochs};
        testingData = data.decoding2.epochs;
    case 3
        epochsForTrain = {data.training1.epochs, data.training2.epochs, data.training3.epochs, ...
                           data.decoding1.epochs, data.decoding2.epochs};
        testingData = data.decoding3.epochs;
    otherwise
        error('Session number must be 1, 2, or 3.');
end

trainingData = combineEpochs(epochsForTrain);

[decoder, modeloutput] = computeDecoderNew(trainingData.data, trainingData.labels, cfg);

%%
[testingData.posteriors, classoutput] = singleClassificationNew(decoder, testingData.data, testingData.labels, 0, ...
                                              decoder.leftElectrodeIndices, decoder.rightElectrodeIndices);


testingData.newlabels = testingData.labels;
testingData.newlabels(testingData.labels == 2) = 1;

figure;

% Define details for each subplot in a structure array:
% For each subplot, specify the label, which field to use ('labels' or 'newlabels'),
% and the desired y-axis limits.
subplotInfo = { ...
    struct('label', 1, 'dataField', 'labels',    'ylim', [0 0.3]), ...
    struct('label', 2, 'dataField', 'labels',    'ylim', [0 0.3]), ...
    struct('label', 0, 'dataField', 'labels',    'ylim', [0 0.3]), ...
    struct('label', 1, 'dataField', 'newlabels', 'ylim', [0 0.5])};

for i = 1:length(subplotInfo)
    subplot(4,1,i);
    cur = subplotInfo{i};
    % Logical index for current label
    idx = testingData.(cur.dataField) == cur.label;
    histogram(testingData.posteriors(idx), ...
        'Normalization', 'probability', ...
        'FaceAlpha', 0.5, ...
        'FaceColor', 'r', ...
        'EdgeColor', 'none', ...
        'NumBins', 10);
    title(sprintf('Label = %d', cur.label));
    xlabel('Posterior Probability');
    ylabel('Count (normalized)');
    xlim([0 1]);
    ylim(cur.ylim);
end

% Add a descriptive overall title to the figure
suptitle(sprintf('Subject %s - Session %d Posterior Probabilities', subjectID, session));

%% Save the Figure
saveFileName = sprintf('%s_session%d.png', subjectID, session);
savePath = fullfile(figpath, saveFileName);
saveas(gcf, savePath);



%%
decoder.decision_threshold = 0.35;
decoder.eegChannels = cfg.eegChannels; 
decoder.eogChannels = cfg.eogChannels;
%decoder.eogFilter = cfg.eogFilter;
decoder.spectralFilter = cfg.spectralFilter;
decoder.threshold = 0.35;
decoder.chantoremove = chanIndices;

decoder.resample.time = decoder.resample.time - decoder.epochOnset;
decoder.psd.time = decoder.psd.time - decoder.epochOnset;
decoder.riemann.time = decoder.riemann.time - decoder.epochOnset;
decoder.asynchronous.threshold = 0.8;
decoder.asynchronous.filterLength = 1;

decoder.cmCV = cmCV;
decoder.accCV = accuracy;
decoder.tnrCV = tnr;
decoder.tprCV = tpr;
decoder.cm = cm;
decoder.NDside = [];
decoder.onlinePosteriors = [];
decoder.subjectID = subjectID;
decoder.datetime = datetime;
disp(' ');
disp('Decoder Updated at');
disp(decoder.datetime);

save(sprintf('./decoders/%s_decoder.mat', subjectID), 'decoder');
save('../cnbiLoop/decoder.mat', 'decoder');
%end
