%function computeModel(subjectID)

%%%%%%%%%%%%%%%%%%%%
%% Initialization %%
%%%%%%%%%%%%%%%%%%%%
clearvars -except subjectID plot_flag;
close all; clc; rng('default');
addpath(genpath('../functions'));
%%%%%%%%%%%%%%%%%%
%% Load dataset %%
%%%%%%%%%%%%%%%%%%
dataPath = [pwd '/../../data/'];
dataInfo = dir([dataPath '/' subjectID '_20*']);

disp(['Loading the data from ' dataInfo.name]);
[calibration, training, decoding] = loadData([dataInfo.folder '/' dataInfo.name '/*']);
%%
cfg = setParams(training.header);

% calibration.trigger = calibration.data(:, cfg.triggerChannel);
% calibration.data = calibration.data(:, [cfg.eegChannels, cfg.eogChannels]);

training.trigger = training.data(:, cfg.triggerChannel);
training.data = training.data(:, [cfg.eegChannels, cfg.eogChannels]);
training.index = computeIndex(training.trigger,4);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Preprocess EEG signals %%                   
%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[b,a] = butter(cfg.spectralFilter.order, cfg.spectralFilter.freqs./(cfg.fsamp/2), 'bandpass');
training.data(:, [cfg.eegChannels, cfg.eogChannels]) = filter(b, a, training.data(:, [cfg.eegChannels, cfg.eogChannels]));


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Plot raw signal %%                   
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% gap = 50;
% figure('Name', 'eog calibration');
% hold on
% subplot(1,2,1);
% %plot(bsxfun(@plus, training.data(1e4:training.eof(1)-2e4, [cfg.eegChannels cfg.eogChannels]), 0:gap:gap*(length([cfg.eegChannels cfg.eogChannels])-1)));
% plot(bsxfun(@plus, training.data(1e4:training.eof(1)-2e4, [cfg.eegChannels]), 0:gap:gap*(length([cfg.eegChannels])-1)));
% yLim = get(gca, 'ylim');
% subplot(1,2,2);
% ex_data = training.data(:, cfg.eegChannels);
% plot(bsxfun(@plus, ex_data(1e4:training.eof-1e4, :), 0:gap:gap*(length(cfg.eegChannels)-1)));
% ylim(yLim);
% legend('Fz', 'Cz', 'Oz')
% xlabel('Time (Samples)')
% ylabel('Amplitude (\muV)')
% title('Raw signal')
% cfg.spectralFilter.b = b;
% cfg.spectralFilter.a = a;
% ylim([0 4000])

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Interpolate noisy channels %%                   
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

num_channels = length(cfg.chanlocs);

% Initialize all channels as good (no interpolation needed)
bad_elec = false(1, num_channels);

%Specify channel name(s) of noisy channel(s)
channels_to_interp = {'AF4'};
bad_channel_index = find(ismember({cfg.chanlocs.labels}, channels_to_interp));
bad_elec(bad_channel_index) = true;


[eeg_interp_data, interpFun] = eeg_interp(training.data, cfg.chanlocs, bad_elec);

% Replace the original data with the interpolated data
training.data = eeg_interp_data;

% %plot
% gap = 50;
% figure('Name', 'eog calibration');
% hold on
% plot(bsxfun(@plus, training.data(1e4:training.eof(1)-2e4, [cfg.eegChannels]), 0:gap:gap*(length([cfg.eegChannels])-1)));
% yLim = get(gca, 'ylim');
% xlabel('Time (Samples)')
% ylabel('Amplitude (\muV)')
% title('Raw signal')



%%%%%%%%%%%%%%
%% Epoching %%
%%%%%%%%%%%%%%
% Define the electrode names for left and right
LeftElectrodes = {'P1', 'P3', 'P5', 'P7', 'PO3', 'PO7'};
RightElectrodes = {'P2', 'P4', 'P6', 'P8', 'PO4', 'PO8'};

% Map electrode names to indices
leftElectrodeIndices = find(ismember({cfg.chanlocs.labels},LeftElectrodes));
rightElectrodeIndices = find(ismember({cfg.chanlocs.labels},RightElectrodes));


% Initialize epochs.data with the correct dimensions
n_samples = length(cfg.epochSamples);
n_channels = length(LeftElectrodes); % Should be 6
n_trials = length(training.index.pos);
epochs.data = nan(n_samples, n_channels, n_trials); % time samples x channels x trials
epochs.labels = training.index.typ;
epochs.file_id = nan(n_trials, 1);
epochs.file_type = cell(n_trials, 1);

for i_trial = 1:n_trials
    label = training.index.typ(i_trial);

    % Select electrode indices based on the label
    if label == 1 || label == 0 % Distractor on right
        electrodeIndices = leftElectrodeIndices;
    elseif label == 2 %|| label == 0 % Distractor on left or no distractor
        electrodeIndices = rightElectrodeIndices;
    else
        error('Unknown label');
    end

    % Extract the data for the selected electrodes
    epochData = training.data(training.index.pos(i_trial) + cfg.epochSamples, electrodeIndices);

    % Store the data in epochs.data
    epochs.data(:, :, i_trial) = epochData;

    % Find the file ID for the current trial
    temp = find(training.index.pos(i_trial) <= training.eof, 1, 'first');
    epochs.file_id(i_trial) = temp;
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Remove bad epochs based on amp %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% rm_index = squeeze(any(any(abs(epochs.data) > 100)));
% epochs.data(:, :, rm_index) = []; 
% epochs.file_id(rm_index) = []; 
% epochs.file_type(rm_index) = [];
% 
% num_bad_epochs = sum(rm_index);
% disp(num_bad_epochs);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Baseline correction %%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
baseline_period = [-0.2, 0];
baseline_indices = find(cfg.epochTime >= baseline_period(1) & cfg.epochTime <= baseline_period(2));
baseline = mean(epochs.data(baseline_indices, :, :), 1); % [1 x channels x trials]
epochs.data = epochs.data - baseline;

%%%%%%%%%%%%%%%%%
%% Calc averages %%
%%%%%%%%%%%%%%%%%

PO7_data_dright = mean(epochs.data(:, 6, epochs.labels == 1), 3);
PO8_data_dleft = mean(epochs.data(:, 6, epochs.labels == 2), 3);
PO8_data_dnone = mean(epochs.data(:, 6, epochs.labels == 0), 3);

grand_average_contra = mean([PO8_data_dleft, PO7_data_dright], 2);
grand_average_dnone = PO8_data_dnone;


%% Plot figure %%
%%%%%%%%%%%%%%%%%

cfg.plotColor = {
    [1, 0, 0],        % Red
    [0, 0, 1],        % Blue
    [1, 0.6, 0.6],    % Light Red
    [0.6, 0.6, 1],    % Light Blue
    [0, 1, 0]         % Green
};

%plot ERP distractor left, right, none
figure('Name','ERP Distractor Lefts vs Right vs No Distractor')
plot(cfg.epochTime, PO7_data_dright, 'Color', cfg.plotColor{1}, cfg.plotOption{:}); hold on;
plot(cfg.epochTime, PO8_data_dleft, 'Color', cfg.plotColor{2}, cfg.plotOption{:}); hold on;
plot(cfg.epochTime, PO8_data_dnone, 'Color', cfg.plotColor{5}, cfg.plotOption{:}); 
xlabel('Time (s)')
ylabel('Amplitude (\muV)');
xline(0);
yline(0);
xlim([-0.1 0.6]); 
ylim([-15 15])
legend('right contra','left contra','none');

%plot ERP distractor left+right vs none
figure('Name','ERP Distractor vs No Distractor')
plot(cfg.epochTime, grand_average_contra, 'Color', cfg.plotColor{1}, cfg.plotOption{:}); hold on;
plot(cfg.epochTime, grand_average_dnone, 'Color', cfg.plotColor{5}, cfg.plotOption{:}); 
xlabel('Time (s)')
ylabel('Amplitude (\muV)');
xline(0);
yline(0);
xlim([-0.1 0.6]); 
ylim([-15 15])
legend('contra', 'none');

% single trial plots
%relabel epochs to be 1 (Distractor) and 0 (No Distractor)
epochs.labels(epochs.labels == 2) = 1;

figure('Name', 'single-trial distractor either left or right vs none')
subplot(2, 1, 1);
imagesc(cfg.epochTime, 1:sum(epochs.labels == 0), squeeze(epochs.data(:, 6, epochs.labels==0))'); caxis([-20 20]);
title('No Distractor')
xlim([-0.1 0.4])
subplot(2, 1, 2);
imagesc(cfg.epochTime, 1:sum(epochs.labels == 1), squeeze(epochs.data(:, 6, epochs.labels == 1))');  caxis([-20 20]);
title('Distractor')
xlim([-0.1 0.4])


figure('Name', 'r-squared')
r2 = compute_r2(epochs.data, epochs.labels);
imagesc(cfg.epochTime, 1:size(epochs.data, 2), r2');
xlim([-0.1 0.7])

drawnow;

%%%%%%%%%%%%%%%%%%%%
%% Classification %%
%%%%%%%%%%%%%%%%%%%%


n_files = length(training.eof);
epochs.posteriors = nan(length(epochs.labels), 1);
for i_file = 1:n_files
    train_index = epochs.file_id ~= i_file; 
    test_index = epochs.file_id == i_file;
    decoder = computeDecoder(epochs.data(:, :, train_index), epochs.labels(train_index), cfg);
    epochs.posteriors(test_index) = singleClassification(decoder, epochs.data(:, :, test_index));
end

disp('== Synchronous Classification == ');
[x, y, t, auc, opt] = perfcurve(~epochs.labels,1-epochs.posteriors, 1, 'Prior', 'uniform');
threshold = 1-t(x == opt(1) & y == opt(2));
disp(['AUC score : ' num2str(auc, '%.2f') ' Threshold: ' num2str(threshold, '%.2f')]);
disp('Confusion Matrix: ');
cm = confusionmat(logical(epochs.labels), (epochs.posteriors >= threshold));
disp(cm);
tnr = cm(1,1) / sum(cm(1, :));
tpr = cm(2,2) / sum(cm(2, :));
disp(['TNR: ' num2str(tnr, '%.2f') ' TPR: ' num2str(tpr, '%.2f')]);

decoder = computeDecoder(epochs.data, epochs.labels, cfg);
epochs.posteriors = singleClassification(decoder, epochs.data);

decoder.decision_threshold = threshold;
decoder.eegChannels = cfg.eegChannels; 
decoder.eogChannels = cfg.eogChannels;
%decoder.eogFilter = cfg.eogFilter;
decoder.spectralFilter = cfg.spectralFilter;
decoder.threshold = threshold;

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