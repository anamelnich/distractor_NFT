subjectID = 'e13';

%%%%%%%%%%%%%%%%%%%%
%% Initialization %%
%%%%%%%%%%%%%%%%%%%%
clearvars -except subjectID  sessions plot_flag;
close all; clc; rng('default');
addpath(genpath('../functions'));
figpath = '../../Figures/';
%%%%%%%%%%%%%%%%%%%%%%
%% Load EEG dataset %%
%%%%%%%%%%%%%%%%%%%%%%
dataPath = [pwd '/../../data/'];
data = loadData(dataPath, subjectID);
delete sopen.mat

% cfg = setParams(data.training1.header);
cfg = setParams(data.eogcalibration1.header);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Load behavioral dataset %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
behData = loadBehData(dataPath, subjectID);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Behavioral Analysis - Distractor Task%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

behfields = fieldnames(behData);
for i = 1:numel(behfields)
    fname = behfields{i};

    if contains(fname,'stroop')
        
        loc = behData.(fname);
        rowsToKeep = ~isnan(loc.Reaction_Time);
        fieldsStroop = fieldnames(loc);
        for k = 1:numel(fieldsStroop)
            fieldStroop = fieldsStroop{k}
            loc.(fieldStroop) = loc.(fieldStroop)(rowsToKeep);
        end
        loc.RTn = loc.Reaction_Time(loc.trial_type==0); %neutral
        loc.RTn_mean = mean(loc.RTn);
        loc.RTn_std = std(loc.RTn);
        loc.RTc = loc.Reaction_Time(loc.trial_type==1); %congruent
        loc.RTc_mean = mean(loc.RTc);
        loc.RTc_std = std(loc.RTc);
        loc.RTi = loc.Reaction_Time(loc.trial_type==2); %incongruent
        loc.RTi_mean = mean(loc.RTi);
        loc.RTi_std = std(loc.RTi);
    else
        loc = behData.(fname);
        loc.RTd = loc.RT(loc.trial_type==1);
        loc.RTd_mean = mean(loc.RTd);
        loc.RTd_std = std(loc.RTd);
        
        loc.RTd_correct = loc.RT(loc.trial_type==1 & loc.response == 1);
        loc.RTd_mean_corr = mean(loc.RTd_correct);
        loc.RTd_std_corr = std(loc.RTd_correct);
        
        loc.RTnd = loc.RT(loc.trial_type==0);
        loc.RTnd_mean = mean(loc.RTnd);
        loc.RTnd_std = std(loc.RTnd);
        
        loc.RTnd_correct = loc.RT(loc.trial_type==0 & loc.response == 1);
        loc.RTnd_mean_corr = mean(loc.RTnd_correct);
        loc.RTnd_std_corr = std(loc.RTnd_correct);
    end 
    behData.(fname) = loc;

end
%%
behfields = fieldnames(behData);
for i = 1:numel(behfields)
    fname = behfields{i};
    if contains(fname,'stroop')
        plotStroopBeh(behData.(fname),[subjectID ' ' fname],figpath);
    elseif contains(fname,'decoding')
    else
        plotBehavior(behData.(fname),cfg,[subjectID ' ' fname],figpath);
    end
end


%%%%%%%%%%%%%%%%%%%%%%%%
%% Set data structure %%                   
%%%%%%%%%%%%%%%%%%%%%%%%

fields = fieldnames(data);
for i = 1:numel(fields)
    fname = fields{i};
    if isempty(data.(fname))
        data = rmfield(data, fname);
        continue;
    end
    data.(fname) = preprocessDataset(data.(fname), cfg,fname);
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Remove non-EEG channels %%                   
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
chanRemove = {'M1','M2','EOG'};
chanIndices = find(ismember(cfg.chanLabels,chanRemove));
cfg.chanLabels(chanIndices) = [];
fields = fieldnames(data);
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
    data.(fname).data(:,:) = filtfilt(b,a,data.(fname).data(:,:));
end

%EOG filter
[b,a] = butter(cfg.EOG.spectralFilter.order, cfg.EOG.spectralFilter.freqs./(cfg.fsamp/2), 'bandpass');
cfg.EOG.spectralFilter.b = b;
cfg.EOG.spectralFilter.a = a;
for i = 1:numel(fields)
    fname = fields{i};
    data.(fname).EOG(:,:) = filtfilt(b,a,data.(fname).EOG(:,:));
end

%%%%%%%%%%%%%%%%%%%%
%% EOG Regression %%                   
%%%%%%%%%%%%%%%%%%%%
% for i = 1:numel(fields)
%     fname = fields{i};
%     if contains(fname, 'eogcalibration')
%         eogFilter = filterEOG(data.(fname).data, data.fname.EOG);
%         cfg.eogFilter = eogFilter;
% 
%     end
% end 
% 
% training.data(:, cfg.eegChannels) = training.data(:, cfg.eegChannels) - training.data(:, cfg.eogChannels) * cfg.eogFilter;
%%%%%%%%%%%%%%%%%%%%%
%% %%%% ICA %%%%%%%%%                   
%%%%%%%%%%%%%%%%%%%%%
for i = 1:numel(fields)
    fname = fields{i};
    dataStruct = data.(fname)
    data.(fname) = processICA(dataStruct);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Remove frontal EEG channels %%                   
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
chanRemove = {'FP1','FP2','FPZ'};
chanIndices = find(ismember(cfg.chanLabels,chanRemove));
cfg.chanLabels(chanIndices) = [];
for i = 1:numel(fields)
    fname = fields{i};
    data.(fname).data(:, chanIndices) = [];
end

%%%%%%%%%%%%%%
%% Epoching %%
%%%%%%%%%%%%%%
for i = 1:numel(fields)
    fname = fields{i};
    dataStruct = data.(fname)
    epochs.data = nan(length(cfg.epochSamples), length(cfg.chanlocs), length(dataStruct.index.pos));
    epochs.EOG = nan(length(cfg.epochSamples), length(cfg.eogChannels), length(dataStruct.index.pos));
    epochs.labels = dataStruct.index.typ;
    epochs.file_id = nan(length(dataStruct.index.typ), 1);
    epochs.file_type = cell(length(dataStruct.index.typ), 1);

    for i_trial = 1:length(dataStruct.index.pos)
        epochs.data(:, :, i_trial) = dataStruct.data(dataStruct.index.pos(i_trial) + cfg.epochSamples, :); 
        epochs.EOG(:, :, i_trial) = dataStruct.EOG(dataStruct.index.pos(i_trial) + cfg.epochSamples, :);
        temp = find(dataStruct.index.pos(i_trial) <= dataStruct.eof, 1, 'first');
        epochs.file_id(i_trial) = temp;
    end
    data.(fname).epochs = epochs;
end 


%%%%%%%%%%%%%%%
%% ERP plots %%
%%%%%%%%%%%%%%%

for i = 1:numel(fields)
    fname = fields{i};
    if contains(fname,'stroop')

    elseif contains(fname,'eogcalibration')

    else
        data.(subjectID).(fname) = computeGrandAvg(data.(fname).epochs.data, ...
            data.(fname).epochs.labels,data.(fname).epochs.file_id,cfg,[subjectID ' ' fname]);
    end 
end

%%%%%%%%%%%%%%%
%% EOG plots %%
%%%%%%%%%%%%%%%

for i = 1:numel(fields)
    fname = fields{i};
    if contains(fname, 'stroop')
        continue;
    elseif contains(fname, 'eogcalibration')
        data.(fname).EOGmVtoDegree = convertEOGtoAngles(data.(fname).epochs.EOG, data.(fname).epochs.labels, cfg, [subjectID ' ' fname], figpath)
    else
        plotEOG(data.(fname).epochs.EOG, data.(fname).epochs.labels, cfg, [subjectID ' ' fname], figpath)
    end
end

%%%%%%%%%%%%%%%
%% Save Data %%
%%%%%%%%%%%%%%%
save([subjectID '.mat'], '-struct', 'data', subjectID);


%%%%%%%%%%%%%%%
%% Load Data %%
%%%%%%%%%%%%%%%
data = load([subjectID '.mat']);


%% RT analysis
behtraining.meanRT = mean(behtraining.RT);
behtraining.stdRT = std(behtraining.RT);

behtraining.meanRT_d = mean(behtraining.RT(epochs.labels == 1 | epochs.labels == 2));
behtraining.stdRT_d = std(behtraining.RT(epochs.labels == 1 | epochs.labels == 2));

behtraining.meanRT_nd = mean(behtraining.RT(epochs.labels == 0));
behtraining.stdRT_nd = std(behtraining.RT(epochs.labels == 0));

figure;
histogram(behtraining.RT(epochs.labels == 1),'FaceColor','red'); hold on;
histogram(behtraining.RT(epochs.labels == 2),'FaceColor','blue'); hold on;
histogram(behtraining.RT(epochs.labels == 0),'FaceColor','green');

%remove trails > 3STDEV RT
RTindices = find((behtraining.RT > (behtraining.meanRT + behtraining.stdRT*3))|(behtraining.RT < (behtraining.meanRT - behtraining.stdRT*3)));
fprintf('n of trials with RTs out of range: %d\n', length(RTindices));
epochs.data(:,:,RTindices)=[];
fieldNames = fieldnames(epochs);
for i = 2:length(fieldNames)
    fieldName = fieldNames{i};
    epochs.(fieldName)(RTindices, :) = [];
end


 %% Averaging for statistics
% timeIdx = 358:461; %0.2-0.6 msec
% 
% avg_offline_dist = mean( offline.ROIerp.grand_average_contra(timeIdx) );
% avg_offline_nodist = mean( offline.ROIerp.grand_average_dnone(timeIdx) );
% 
% avg_online1_dist = mean( online1.ROIerp.grand_average_contra(timeIdx) );
% avg_online1_nodist = mean( online1.ROIerp.grand_average_dnone(timeIdx) );
% 
% avg_online2_dist = mean( online2.ROIerp.grand_average_contra(timeIdx) );
% avg_online2_nodist = mean( online2.ROIerp.grand_average_dnone(timeIdx) );
% 
% subjectRow = [ ...
%     avg_offline_nodist, ...
%     avg_offline_dist, ...
%     avg_online1_nodist, ...
%     avg_online1_dist, ...
%     avg_online2_nodist, ...
%     avg_online2_dist ];
% load('ERPamplitudes.mat');
% ERPamplitudes(8, :) = subjectRow;
% save('ERPamplitudes.mat', 'ERPamplitudes');



%%
topoplotxtime(offline.allElect.allElect_dright, offline.allElect.allElect_dleft, offline.allElect.allElect_nd, cfg.chanlocs, ...
    0.6, 1.1, 0.1, 0.05,'Offline') %roughly 100 to 0.6
%%
topoplotxtime(online1.allElect.allElect_dright, online1.allElect.allElect_dleft, online1.allElect.allElect_nd, cfg.chanlocs, ...
    0.6, 1.1, 0.1, 0.05,'Online 1')
%%
topoplotxtime(online2.allElect.allElect_dright, online2.allElect.allElect_dleft, online2.allElect.allElect_nd, cfg.chanlocs, ...
    0.6, 1.1, 0.1, 0.05, 'Online 2')
%%
topoplotxtime(grandaverage_dright_topo.offline, grandaverage_dleft_topo.offline, ...
    grandaverage_dnone_topo.offline, cfg.chanlocs, ...
    0.6, 1.1, 0.1, 0.05,'Offline Grand Average')
%%
topoplotxtime(grandaverage_dright_topo.online1, grandaverage_dleft_topo.online1, ...
    grandaverage_dnone_topo.online1, cfg.chanlocs, ...
    0.6, 1.1, 0.1, 0.05,'Online 1 Grand Average')
%%
topoplotxtime(grandaverage_dright_topo.online2, grandaverage_dleft_topo.online2, ...
    grandaverage_dnone_topo.online2, cfg.chanlocs, ...
    0.6, 1.1, 0.1, 0.05,'Online 2 Grand Average')


%%
dright_diff = grandaverage_dright_topo.online1 - grandaverage_dright_topo.offline;
dleft_diff = grandaverage_dleft_topo.online1 - grandaverage_dleft_topo.offline;
dnone_diff = grandaverage_dnone_topo.online1 - grandaverage_dnone_topo.offline;
%%
topoplotxtime(dright_diff, dleft_diff, ...
    dnone_diff, cfg.chanlocs, ...
    0.6, 1.1, 0.1, 0.05,'Grand Average Difference (Online1 - Offline)')


%% single trial plots
%relabel epochs to be 1 (Distractor) and 0 (No Distractor)
%epochs.labels(epochs.labels == 2) = 1;
% 
% figure('Name', 'single-trial distractor either left or right vs none')
% subplot(2, 1, 1);
% imagesc(cfg.epochTime, 1:sum(epochs.labels == 0), squeeze(epochs.data(:, 6, epochs.labels==0))'); caxis([-20 20]);
% title('No Distractor')
% xlim([-0.1 0.4])
% subplot(2, 1, 2);
% imagesc(cfg.epochTime, 1:sum(epochs.labels == 1), squeeze(epochs.data(:, 6, epochs.labels == 1))');  caxis([-20 20]);
% title('Distractor')
% xlim([-0.1 0.4])
% 
% 
% figure('Name', 'r-squared')
% r2 = compute_r2(epochs.data, epochs.labels);
% imagesc(cfg.epochTime, 1:size(epochs.data, 2), r2');
% xlim([-0.1 0.7])
% 
% drawnow;

%%%%%%%%%%%%%%%%%%%%
%% Classification %%
%%%%%%%%%%%%%%%%%%%%


n_files = length(training.eof);
offline.posteriors = nan(length(offline.labels), 1);
for i_file = 1:n_files
    train_index = offline.file_id ~= i_file; 
    test_index = offline.file_id == i_file;
    decoder = computeDecoder(offline.preprocEpoch(:, :, train_index), offline.labels(train_index), cfg);
    offline.posteriors(test_index) = singleClassification(decoder, offline.preprocEpoch(:, :, test_index), offline.labels(test_index),0,decoder.leftElectrodeIndices,decoder.rightElectrodeIndices);
end
offline.newlabels = offline.labels;
offline.newlabels(offline.labels == 2) = 1;
disp('== Synchronous Classification == ');
% [x, y, t, auc, opt] = perfcurve(~epochs.labels,1-epochs.posteriors, 1, 'Prior', 'uniform');
% threshold = 1-t(x == opt(1) & y == opt(2));
[x, y, t, auc, opt] = perfcurve(offline.newlabels,offline.posteriors, 1, 'Prior', 'uniform');
% threshold = t(x == opt(1) & y == opt(2));
threshold = 0.5
disp(['AUC score : ' num2str(auc, '%.2f') ' Threshold: ' num2str(threshold, '%.2f')]);
disp('Confusion Matrix: ');
cm = confusionmat(logical(offline.newlabels), (offline.posteriors >= threshold));
disp(cm);
tnr = cm(1,1) / sum(cm(1, :));
tpr = cm(2,2) / sum(cm(2, :));
disp(['TNR: ' num2str(tnr, '%.2f') ' TPR: ' num2str(tpr, '%.2f')]);
TN_t = cm(1, 1); 
TP_t = cm(2, 2); 
cmTotal_t = sum(cm(:));
cmCorrect_t = TP_t + TN_t;
accuracy_t = (cmCorrect_t / cmTotal_t) * 100;
disp(accuracy_t);
%%
% decoder1 = computeDecoder(decoding1_epochs.data, decoding1_epochs.labels, cfg);
decoding1_epochs.threshold = threshold;
decoding1_epochs.posteriors = singleClassification(decoder, decoding1_epochs.data, decoding1_epochs.labels, 0, decoder.leftElectrodeIndices,decoder.rightElectrodeIndices);
decoding1_epochs.newlabels = decoding1_epochs.labels;
decoding1_epochs.newlabels(decoding1_epochs.labels==2)=1;
disp(['AUC score : ' num2str(auc, '%.2f') ' Threshold: ' num2str(threshold, '%.2f')]);
disp('Confusion Matrix: ');
cm = confusionmat(logical(decoding1_epochs.newlabels), (decoding1_epochs.posteriors >= threshold));
disp(cm);
tnr = cm(1,1) / sum(cm(1, :));
tpr = cm(2,2) / sum(cm(2, :));
disp(['TNR: ' num2str(tnr, '%.2f') ' TPR: ' num2str(tpr, '%.2f')]);
TN_t = cm(1, 1); 
TP_t = cm(2, 2); 
cmTotal_t = sum(cm(:));
cmCorrect_t = TP_t + TN_t;
accuracy_t = (cmCorrect_t / cmTotal_t) * 100;
disp(accuracy_t);
%%
decoding2_epochs.threshold = threshold;
decoding2_epochs.posteriors = singleClassification(decoder, decoding2_epochs.data, decoding2_epochs.labels, 0, decoder.leftElectrodeIndices,decoder.rightElectrodeIndices);
decoding2_epochs.newlabels = decoding2_epochs.labels;
decoding2_epochs.newlabels(decoding2_epochs.labels==2)=1;
disp(['AUC score : ' num2str(auc, '%.2f') ' Threshold: ' num2str(threshold, '%.2f')]);
disp('Confusion Matrix: ');
cm = confusionmat(logical(decoding2_epochs.newlabels), (decoding2_epochs.posteriors >= threshold));
disp(cm);
tnr = cm(1,1) / sum(cm(1, :));
tpr = cm(2,2) / sum(cm(2, :));
disp(['TNR: ' num2str(tnr, '%.2f') ' TPR: ' num2str(tpr, '%.2f')]);
TN_t = cm(1, 1); 
TP_t = cm(2, 2); 
cmTotal_t = sum(cm(:));
cmCorrect_t = TP_t + TN_t;
accuracy_t = (cmCorrect_t / cmTotal_t) * 100;
disp(accuracy_t);

%% Build model based on both offline and online 1

newepochs = struct();
newepochs.data = cat(3, epochs1.data, depochs1.data);
newepochs.labels = [epochs1.labels; depochs1.labels];
newepochs.file_id = [epochs1.file_id; depochs1.file_id];

n_files = length(training.eof);
newepochs.posteriors = nan(length(newepochs.labels), 1);
for i_file = 1:n_files
    train_index = newepochs.file_id ~= i_file; 
    test_index = newepochs.file_id == i_file;
    decoder = computeDecoder(newepochs.data(:, :, train_index), newepochs.labels(train_index), cfg);
    newepochs.posteriors(test_index) = singleClassification(decoder, newepochs.data(:, :, test_index), newepochs.labels(test_index),0,decoder.leftElectrodeIndices,decoder.rightElectrodeIndices);
end
newepochs.newlabels = newepochs.labels;
newepochs.newlabels(newepochs.labels == 2) = 1;
disp('== Synchronous Classification == ');
% [x, y, t, auc, opt] = perfcurve(~epochs.labels,1-epochs.posteriors, 1, 'Prior', 'uniform');
% threshold = 1-t(x == opt(1) & y == opt(2));
[x, y, t, auc, opt] = perfcurve(newepochs.newlabels,newepochs.posteriors, 1, 'Prior', 'uniform');
% threshold = t(x == opt(1) & y == opt(2));
threshold = 0.5
disp(['AUC score : ' num2str(auc, '%.2f') ' Threshold: ' num2str(threshold, '%.2f')]);
disp('Confusion Matrix: ');
cm = confusionmat(logical(newepochs.newlabels), (newepochs.posteriors >= threshold));
disp(cm);
tnr = cm(1,1) / sum(cm(1, :));
tpr = cm(2,2) / sum(cm(2, :));
disp(['TNR: ' num2str(tnr, '%.2f') ' TPR: ' num2str(tpr, '%.2f')]);
%%
decoder2 = computeDecoder(newepochs.data, newepochs.labels, cfg);
decoder2.threshold = threshold;
depochs2.posteriors = singleClassification(decoder2, depochs2.data, depochs2.labels, 0, decoder.leftElectrodeIndices,decoder.rightElectrodeIndices);
% epochs2.posteriors = singleClassification(decoder2, epochs2.data, epochs2.labels, 0, decoder.leftElectrodeIndices,decoder.rightElectrodeIndices);


%% Plot class hist
figure('Name', 'Posterior Probabilities Distribution');
histogram(epochs1.posteriors(epochs1.labels ==0), 'BinWidth', 0.05, 'FaceColor', 'green', 'FaceAlpha', 0.5); hold on;
histogram(epochs1.posteriors(epochs1.labels ==1), 'BinWidth', 0.05, 'FaceColor', 'red', 'FaceAlpha', 0.5); hold on;
histogram(epochs1.posteriors(epochs1.labels ==2), 'BinWidth', 0.05, 'FaceColor', 'blue', 'FaceAlpha', 0.5); 
xline(threshold, 'k--', 'LineWidth', 2);
xlabel('Posterior Probability');
ylabel('Frequency');
title('Offline Pre');
legend('No Distractor', 'Distractor Right','Distractor Left');
ylim([0 20]);
set(gca, 'FontSize', 16);

figure('Name', 'Posterior Probabilities Distribution');
histogram(depochs1.posteriors(depochs1.labels ==0), 'BinWidth', 0.05, 'FaceColor', 'green', 'FaceAlpha', 0.5); hold on;
histogram(depochs1.posteriors(depochs1.labels ==1), 'BinWidth', 0.05, 'FaceColor', 'red', 'FaceAlpha', 0.5); hold on;
histogram(depochs1.posteriors(depochs1.labels ==2), 'BinWidth', 0.05, 'FaceColor', 'blue', 'FaceAlpha', 0.5); 
xline(threshold, 'k--', 'LineWidth', 2);
xlabel('Posterior Probability');
ylabel('Frequency');
title('Online 1');
legend('No Distractor', 'Distractor Right','Distractor Left');
ylim([0 20]);
set(gca, 'FontSize',16);

figure('Name', 'Posterior Probabilities Distribution');
histogram(depochs2.posteriors(depochs2.labels ==0), 'BinWidth', 0.05, 'FaceColor', 'green', 'FaceAlpha', 0.5); hold on;
histogram(depochs2.posteriors(depochs2.labels ==1), 'BinWidth', 0.05, 'FaceColor', 'red', 'FaceAlpha', 0.5); hold on;
histogram(depochs2.posteriors(depochs2.labels ==2), 'BinWidth', 0.05, 'FaceColor', 'blue', 'FaceAlpha', 0.5); 
xline(threshold, 'k--', 'LineWidth', 2);
xlabel('Posterior Probability');
ylabel('Frequency');
title('Online 2');
legend('No Distractor', 'Distractor Right','Distractor Left');
ylim([0 30]);
set(gca, 'FontSize',16);

figure('Name', 'Posterior Probabilities Distribution');
histogram(epochs2.posteriors(epochs2.labels ==0), 'BinWidth', 0.05, 'FaceColor', 'green', 'FaceAlpha', 0.5); hold on;
histogram(epochs2.posteriors(epochs2.labels ==1), 'BinWidth', 0.05, 'FaceColor', 'red', 'FaceAlpha', 0.5); hold on;
histogram(epochs2.posteriors(epochs2.labels ==2), 'BinWidth', 0.05, 'FaceColor', 'blue', 'FaceAlpha', 0.5); 
xline(threshold, 'k--', 'LineWidth', 2);
xlabel('Posterior Probability');
ylabel('Frequency');
title('Offline Post');
legend('No Distractor', 'Distractor Right','Distractor Left');
ylim([0 15]);
set(gca, 'FontSize',16);

%% class hist separated
figure('Name', 'Posterior Probabilities Distribution');

% First Subplot: No Distractor
subplot(3,1,1);  % Divide the figure into a 3x1 grid and use the first section
histogram(epochs.posteriors(epochs.labels == 0), 'BinWidth', 0.05, ...
    'FaceColor', 'green', 'FaceAlpha', 0.5);
xline(threshold, 'k--', 'LineWidth', 2);
ylabel('Frequency');
title('No Distractor');
ylim([0 40]);
xlim([0 1]);

% Second Subplot: Distractor Right
subplot(3,1,2);  % Use the second section of the 3x1 grid
histogram(epochs.posteriors(epochs.labels == 1), 'BinWidth', 0.05, ...
    'FaceColor', 'red', 'FaceAlpha', 0.5);
xline(threshold, 'k--', 'LineWidth', 2);
ylabel('Frequency');
title('Distractor Right');
ylim([0 40]);
xlim([0 1]);

% Third Subplot: Distractor Left
subplot(3,1,3);  % Use the third section of the 3x1 grid
histogram(epochs.posteriors(epochs.labels == 2), 'BinWidth', 0.05, ...
    'FaceColor', 'blue', 'FaceAlpha', 0.5);
xline(threshold, 'k--', 'LineWidth', 2);
xlabel('Posterior Probability');
ylabel('Frequency');
title('Distractor Left');
ylim([0 40]);
xlim([0 1]);


%% Plot confusion matrix 
cmChart_t = confusionchart(logical(epochs1.labels), (epochs1.posteriors >= decoder1.threshold), ...
    'Normalization', 'row-normalized', 'Title', 'Confusion Matrix Training','FontSize',20); 
%%
depochs1.newlabels = depochs1.labels;
depochs1.newlabels(depochs1.labels == 2) = 1;
cmChart_d1 = confusionchart(logical(depochs1.newlabels), (depochs1.posteriors >= decoder1.threshold), ...
    'Normalization', 'row-normalized', 'Title', 'Confusion Matrix Decoding 1','FontSize',20);
%%
depochs2.newlabels = depochs2.labels;
depochs2.newlabels(depochs2.labels == 2) = 1;
cmChart_d2 = confusionchart(logical(depochs2.newlabels), (depochs2.posteriors >= decoder2.threshold), ...
    'Normalization', 'row-normalized', 'Title', 'Confusion Matrix Decoding 2','FontSize',20);
%% Accuracy

cm_t = confusionmat(logical(epochs1.newlabels), (epochs1.posteriors >= decoder1.threshold));
TN_t = cm_t(1, 1); 
FP_t = cm_t(1, 2); 
FN_t = cm_t(2, 1); 
TP_t = cm_t(2, 2); 
cmTotal_t = sum(cm_t(:));
cmCorrect_t = TP_t + TN_t;
accuracy_t = (cmCorrect_t / cmTotal_t) * 100;
precision_t = TP_t / (TP_t + FP_t);
sensitivity_t = TP_t / (TP_t + FN_t); % Ability to correctly identify positive cases
specificity_t = TN_t / (TN_t + FP_t); % Ability to correctly identify negative cases

cm_d = confusionmat(logical(depochs1.labels), (depochs1.posteriors >= decoder1.threshold));
TN_d = cm_d(1, 1); 
FP_d = cm_d(1, 2); 
FN_d = cm_d(2, 1); 
TP_d = cm_d(2, 2); 
cmTotal_d = sum(cm_d(:));
cmCorrect_d = TP_d + TN_d;
accuracy_d = (cmCorrect_d / cmTotal_d) * 100;
precision_d = TP_d / (TP_d + FP_d);
sensitivity_d = TP_d / (TP_d + FN_d); % Ability to correctly identify positive cases
specificity_d = TN_d / (TN_d + FP_d); % Ability to correctly identify negative cases
% 
cm_d2 = confusionmat(logical(depochs2.newlabels), (depochs2.posteriors >= decoder2.threshold));
TN_d2 = cm_d2(1, 1); 
FP_d2 = cm_d2(1, 2); 
FN_d2 = cm_d2(2, 1); 
TP_d2 = cm_d2(2, 2); 
cmTotal_d2 = sum(cm_d2(:));
cmCorrect_d2 = TP_d2 + TN_d2;
accuracy_d2 = (cmCorrect_d2 / cmTotal_d2) * 100;
precision_d2 = TP_d2 / (TP_d2 + FP_d2);
sensitivity_d2 = TP_d2 / (TP_d2 + FN_d2); % Ability to correctly identify positive cases
specificity_d2 = TN_d2 / (TN_d2 + FP_d2); % Ability to correctly identify negative cases

accuracy = [accuracy_t, accuracy_d, accuracy_d2];
precision = [precision_t, precision_d, precision_d2] * 100; 
sensitivity = [sensitivity_t, sensitivity_d, sensitivity_d2] * 100; 
specificity = [specificity_t, specificity_d, specificity_d2] * 100; 

dataset_labels = {'Offline', 'Online 1', 'Online 2'};
metrics_matrix = [
    accuracy;
    precision;
    sensitivity;
    specificity
];
metric_labels = {'Accuracy', 'Precision', 'Sensitivity', 'Specificity'};

% Create grouped bar chart
figure;
bar(metrics_matrix');
set(gca, 'XTickLabel', dataset_labels, 'FontSize', 12);
ylabel('Percentage (%)', 'FontSize', 14);
legend(metric_labels, 'Location', 'northwest');
title('Performance Metrics Across Datasets', 'FontSize', 16);
ylim([40 90]);

%% Accuracy over runs
uniqueRuns1 = unique(depochs1.file_id);
n_runs1 = length(uniqueRuns1);
accuracy_per_run1 = zeros(n_runs1, 1);

for i = 1:n_runs1
    run_id = uniqueRuns1(i);
    run_indices = (depochs1.file_id == run_id);
    true_labels = logical(depochs1.newlabels(run_indices));
    predicted_labels = depochs1.posteriors(run_indices) >= decoder1.threshold;
    
    cm = confusionmat(true_labels, predicted_labels);
    TN = cm(1, 1);
    FP = cm(1, 2);
    FN = cm(2, 1);
    TP = cm(2, 2);
    
    total = TP + TN + FP + FN;
    correct = TP + TN;
    accuracy_per_run1(i) = (correct / total) * 100;
end

% For depochs2
uniqueRuns2 = unique(depochs2.file_id);
n_runs2 = length(uniqueRuns2);
accuracy_per_run2 = zeros(n_runs2, 1);

for i = 1:n_runs2
    run_id = uniqueRuns2(i);
    run_indices = (depochs2.file_id == run_id);
    true_labels = logical(depochs2.newlabels(run_indices));
    predicted_labels = depochs2.posteriors(run_indices) >= decoder2.threshold;
    
    cm = confusionmat(true_labels, predicted_labels);
    TN = cm(1, 1);
    FP = cm(1, 2);
    FN = cm(2, 1);
    TP = cm(2, 2);
    
    total = TP + TN + FP + FN;
    correct = TP + TN;
    accuracy_per_run2(i) = (correct / total) * 100;
end

% Plot both on the same axes using their actual run IDs
figure; hold on;
plot(uniqueRuns1, accuracy_per_run1, '-o', 'LineWidth', 2, 'MarkerSize', 8, 'DisplayName', 'Online 1 Accuracy');
plot(uniqueRuns2, accuracy_per_run2, '-s', 'LineWidth', 2, 'MarkerSize', 8, 'DisplayName', 'Online 2 Accuracy');

xlabel('Run Number', 'FontSize', 12);
ylabel('Accuracy (%)', 'FontSize', 12);
title('Accuracy Over Time (Per Run)', 'FontSize', 14);
yline(50, 'k--', 'DisplayName', 'Chance Level');
yline(accuracy_t, 'r--', 'DisplayName', 'Offline Accuracy');
legend('Location', 'best');
grid on;
ylim([40 90]);

% Set xticks to cover the entire range (including runs 1-4 and then 5+)
allRuns = [uniqueRuns1; uniqueRuns2];
xticks(min(allRuns):max(allRuns));


%% Save decoder info

decoder.decision_threshold = threshold;
decoder.eegChannels = cfg.eegChannels; 
decoder.eogChannels = cfg.eogChannels;
%decoder.eogFilter = cfg.eogFilter;
decoder.spectralFilter = cfg.spectralFilter;
% decoder.threshold = threshold;

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
