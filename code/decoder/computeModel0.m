function computeModel(subjectID)

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

%%%%%%%%%%%%%%%%%%%%%%%
%% Calculate RT time %%                   
%%%%%%%%%%%%%%%%%%%%%%%

training = computeReactionTimes(training, training.index, cfg.fsamp);
disp(training.RT);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Preprocess EEG signals %%                   
%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[b,a] = butter(cfg.spectralFilter.order, cfg.spectralFilter.freqs./(cfg.fsamp/2), 'bandpass');
training.data(:, [cfg.eegChannels, cfg.eogChannels]) = filter(b, a, training.data(:, [cfg.eegChannels, cfg.eogChannels]));
%eogFilter = filterEOG(calibration.data(calibration.start_index:calibration.end_index, cfg.eegChannels), calibration.data(calibration.start_index:calibration.end_index, cfg.eogChannels));
%cfg.eogFilter = eogFilter;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Plot signals before interp %%                   
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

gap = 50;
figure('Name', 'eog calibration');
hold on
subplot(1,2,1);
%plot(bsxfun(@plus, training.data(1e4:training.eof(1)-2e4, [cfg.eegChannels cfg.eogChannels]), 0:gap:gap*(length([cfg.eegChannels cfg.eogChannels])-1)));
plot(bsxfun(@plus, training.data(1e4:training.eof(1)-2e4, [cfg.eegChannels]), 0:gap:gap*(length([cfg.eegChannels])-1)));
yLim = get(gca, 'ylim');
subplot(1,2,2);
ex_data = training.data(:, cfg.eegChannels);
plot(bsxfun(@plus, ex_data(1e4:training.eof-1e4, :), 0:gap:gap*(length(cfg.eegChannels)-1)));
ylim(yLim);
legend('Fz', 'Cz', 'Oz')
xlabel('Time (Samples)')
ylabel('Amplitude (\muV)')
title('Raw signal')
cfg.spectralFilter.b = b;
cfg.spectralFilter.a = a;
ylim([0 4000])


% gap = 50;
% figure('Name', 'eog calibration');
% subplot(1,2,1);
% plot(bsxfun(@plus, calibration.data(calibration.start_index:calibration.end_index, [parameters.eegChannels parameters.eogChannels]), 0:gap:gap*(length([parameters.eegChannels parameters.eogChannels])-1)));
% yLim = get(gca, 'ylim');
% subplot(1,2,2);
% ex_data = calibration.data(:, parameters.eegChannels) - calibration.data(:, parameters.eogChannels) * eogFilter;
% plot(bsxfun(@plus, ex_data(calibration.start_index:calibration.end_index, :), 0:gap:gap*(length(parameters.eegChannels)-1)));
% ylim(yLim);

% training.data(:, [cfg.eegChannels, cfg.eogChannels]) = filter(b, a, training.data(:, [cfg.eegChannels, cfg.eogChannels]));
%training.data(:, cfg.eegChannels) = training.data(:, cfg.eegChannels) - training.data(:, cfg.eogChannels) * cfg.eogFilter;

% gap = 50;
% figure('Name', 'eog calibration');
% % subplot(1,2,1);
% plot(bsxfun(@plus, training.data(1e4:training.eof(1)-2e4, [parameters.eegChannels parameters.eogChannels]), 0:gap:gap*(length([parameters.eegChannels parameters.eogChannels])-1)));
% yLim = get(gca, 'ylim');
% % subplot(1,2,2);
% % ex_data = training.data(:, cfg.eegChannels) - training.data(:, cfg.eogChannels) * eogFilter;
% % plot(bsxfun(@plus, ex_data(1e4:training.eof-1e4, :), 0:gap:gap*(length(cfg.eegChannels)-1)));
% legend('Fz', 'Cz', 'Oz')
% % ylim(yLim);
% xlabel('Time (Samples)')
% ylabel('Amplitude (\muV)')
% title('Bandpass filtered signal')

% cfg.spectralFilter.b = b;
% cfg.spectralFilter.a = a;

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


gap = 50;
figure('Name', 'eog calibration');
hold on
plot(bsxfun(@plus, training.data(1e4:training.eof(1)-2e4, [cfg.eegChannels]), 0:gap:gap*(length([cfg.eegChannels])-1)));
yLim = get(gca, 'ylim');
xlabel('Time (Samples)')
ylabel('Amplitude (\muV)')
title('Raw signal')



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
epochs.RT = training.RT.typ;
epochs.file_id = nan(n_trials, 1);
epochs.file_type = cell(n_trials, 1);

for i_trial = 1:n_trials
    label = training.index.typ(i_trial);

    % Select electrode indices based on the label
    if label == 1 %|| label == 0  % Distractor on right or no distractor
        electrodeIndices = leftElectrodeIndices;
    elseif label == 2 || label == 0 % Distractor on left or no distractor
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


%%%%%%%%%%%%%%%%%
%% Calc averages %%
%%%%%%%%%%%%%%%%%
PO7 = find(strcmp({cfg.chanlocs.labels}, 'PO7'));
PO8 = find(strcmp({cfg.chanlocs.labels}, 'PO8'));


baseline_period = [-0.2, 0];
baseline_indices = find(cfg.epochTime >= baseline_period(1) & cfg.epochTime <= baseline_period(2));
baseline = mean(epochs.data(baseline_indices, :, :), 1); % [1 x channels x trials]
epochs.data = epochs.data - baseline;
%%
RT = 0.626
PO7_data_dright = mean(epochs.data(:, PO7, epochs.labels == 1), 3);
PO8_data_dright = mean(epochs.data(:, PO8, epochs.labels == 1), 3);
PO7_data_dleft = mean(epochs.data(:, PO7, epochs.labels == 2), 3);
PO8_data_dleft = mean(epochs.data(:, PO8, epochs.labels == 2), 3);
PO7_data_dnone = mean(epochs.data(:, PO7, epochs.labels == 0), 3);
PO8_data_dnone = mean(epochs.data(:, PO8, epochs.labels == 0), 3);

grand_average_ipsi = mean([PO7_data_dleft, PO8_data_dright], 2);
grand_average_contra = mean([PO8_data_dleft, PO7_data_dright], 2);
grand_average_dnone = mean([PO7_data_dnone,PO8_data_dnone],2);

dright_channel_difference = PO7_data_dright - PO8_data_dright;
dleft_channel_difference = PO8_data_dleft - PO7_data_dleft;
dnone_channel_difference = PO7_data_dnone - PO8_data_dnone;
grand_average_channel_difference = mean([dright_channel_difference, dleft_channel_difference], 2);

singletrial_difference_dright = epochs.data(:, PO7, epochs.labels == 1) - epochs.data(:, PO8, epochs.labels == 1);
singletrial_difference_dleft = epochs.data(:, PO8, epochs.labels == 2) - epochs.data(:, PO7, epochs.labels == 2);
singletrial_difference_dnone = epochs.data(:, PO7, epochs.labels == 0) - epochs.data(:, PO8, epochs.labels == 0);

%plot based on RT (below or above average)

PO7_data_dright_lRT = mean(epochs.data(:, PO7, ((epochs.labels == 1) & (epochs.RT == 1))), 3);
PO7_data_dright_hRT = mean(epochs.data(:, PO7, ((epochs.labels == 1) & (epochs.RT == 0))), 3);
PO8_data_dright_lRT = mean(epochs.data(:, PO8, ((epochs.labels == 1) & (epochs.RT == 1))), 3);
PO8_data_dright_hRT = mean(epochs.data(:, PO8, ((epochs.labels == 1) & (epochs.RT == 0))), 3);

PO7_data_dleft_lRT = mean(epochs.data(:, PO7, ((epochs.labels == 2) & (epochs.RT == 1))), 3);
PO7_data_dleft_hRT = mean(epochs.data(:, PO7, ((epochs.labels == 2) & (epochs.RT == 0))), 3);
PO8_data_dleft_lRT = mean(epochs.data(:, PO8, ((epochs.labels == 2) & (epochs.RT == 1))), 3);
PO8_data_dleft_hRT = mean(epochs.data(:, PO8, ((epochs.labels == 2) & (epochs.RT == 0))), 3);

PO7_data_dnone_lRT = mean(epochs.data(:, PO7, ((epochs.labels == 0) & (epochs.RT == 1))), 3);
PO7_data_dnone_hRT = mean(epochs.data(:, PO7, ((epochs.labels == 0) & (epochs.RT == 0))), 3);
PO8_data_dnone_lRT = mean(epochs.data(:, PO8, ((epochs.labels == 0) & (epochs.RT == 1))), 3);
PO8_data_dnone_hRT = mean(epochs.data(:, PO8, ((epochs.labels == 0) & (epochs.RT == 0))), 3);

grand_average_ipsi_lRT = mean([PO7_data_dleft_lRT, PO8_data_dright_lRT], 2);
grand_average_contra_lRT = mean([PO8_data_dleft_lRT, PO7_data_dright_lRT], 2);
grand_average_ipsi_hRT = mean([PO7_data_dleft_hRT, PO8_data_dright_hRT], 2);
grand_average_contra_hRT = mean([PO8_data_dleft_hRT, PO7_data_dright_hRT], 2);
grand_average_dnone_lRT = mean([PO7_data_dnone_lRT, PO8_data_dnone_lRT], 2);
grand_average_dnone_hRT = mean([PO8_data_dleft_hRT, PO7_data_dright_hRT], 2);


%% Plot figure %%
%%%%%%%%%%%%%%%%%
RT=0.585
%plot ERP waveforms

ax = gobjects(3,1)
figure('Name', 'grand-average ERP');
ax(1)=subplot(3, 1, 1);
plot(cfg.epochTime, PO8_data_dleft, 'Color', cfg.plotColor{1}, cfg.plotOption{:}); hold on;
plot(cfg.epochTime, PO7_data_dleft, 'Color', cfg.plotColor{2}, cfg.plotOption{:}); hold on;
title('Distractor on the Left');
ax(2)=subplot(3, 1, 2);
plot(cfg.epochTime, PO7_data_dright, 'Color', cfg.plotColor{1}, cfg.plotOption{:}); hold on;
plot(cfg.epochTime, PO8_data_dright, 'Color', cfg.plotColor{2}, cfg.plotOption{:}); 
title('Distractor on the Right');
ax(3)=subplot(3, 1, 3);
plot(cfg.epochTime, grand_average_contra, 'Color', cfg.plotColor{1}, cfg.plotOption{:}); hold on;
plot(cfg.epochTime, grand_average_ipsi, 'Color', cfg.plotColor{2}, cfg.plotOption{:}); 
title('Average Left and Right');
xlabel('Time (s)')
for i = 1:length(ax)
    axes(ax(i)); 
    ylabel('Amplitude (\muV)');
    xline(0); 
    xline(RT, '--', 'LineWidth', 2); 
    yline(0);
    xlim([-0.1 0.55]); 
    ylim([-10 10])
    legend('contra', 'ipsi');
end

%%
cfg.plotColor = {
    [1, 0, 0],        % Red
    [0, 0, 1],        % Blue
    [1, 0.6, 0.6],    % Light Red
    [0.6, 0.6, 1],    % Light Blue
    [0, 1, 0]         % Green
};
figure('Name', 'grand-average ERP Distractor Right');
plot(cfg.epochTime, PO7_data_dright_lRT, 'Color', cfg.plotColor{1}, cfg.plotOption{:}); hold on; 
plot(cfg.epochTime, PO8_data_dright_lRT, 'Color', cfg.plotColor{2}, cfg.plotOption{:}); hold on; 
plot(cfg.epochTime, PO7_data_dright_hRT, 'Color', cfg.plotColor{3}, cfg.plotOption{:});hold on; 
plot(cfg.epochTime, PO8_data_dright_hRT, 'Color', cfg.plotColor{4}, cfg.plotOption{:});
plot(cfg.epochTime, grand_average_dnone, 'Color', cfg.plotColor{5}, cfg.plotOption{:}); 
xlabel('Time (s)')
ylabel('Amplitude (\muV)');
xline(0);
yline(0);
xlim([-0.1 0.55]); 
ylim([-10 10])
legend('contra low RT', 'ipsi low RT','contra high RT', 'ipsi high RT','no distractor');

figure('Name', 'grand-average ERP Distractor Left');
plot(cfg.epochTime, PO8_data_dleft_lRT, 'Color', cfg.plotColor{1}, cfg.plotOption{:}); hold on; 
plot(cfg.epochTime, PO7_data_dleft_lRT, 'Color', cfg.plotColor{2}, cfg.plotOption{:}); hold on; 
plot(cfg.epochTime, PO8_data_dleft_hRT, 'Color', cfg.plotColor{3}, cfg.plotOption{:}); hold on;
plot(cfg.epochTime, PO7_data_dleft_hRT, 'Color', cfg.plotColor{4}, cfg.plotOption{:}); 
plot(cfg.epochTime, grand_average_dnone, 'Color', cfg.plotColor{5}, cfg.plotOption{:}); 
xlabel('Time (s)')
ylabel('Amplitude (\muV)');
xline(0);
yline(0);
xlim([-0.1 0.55]); 
ylim([-10 10])
legend('contra low RT', 'ipsi low RT','contra high RT', 'ipsi high RT','none');

figure('Name', 'grand-average ERP Distractor Left + Right by RT');
plot(cfg.epochTime, grand_average_contra_lRT, 'Color', cfg.plotColor{1}, cfg.plotOption{:}, 'LineStyle','-'); hold on; 
plot(cfg.epochTime, grand_average_ipsi_lRT, 'Color', cfg.plotColor{2}, cfg.plotOption{:}, 'LineStyle','-'); hold on; 
plot(cfg.epochTime, grand_average_contra_hRT, 'Color', cfg.plotColor{3}, cfg.plotOption{:},'LineStyle','-'); hold on;
plot(cfg.epochTime, grand_average_ipsi_hRT, 'Color', cfg.plotColor{4}, cfg.plotOption{:},'LineStyle','-');
plot(cfg.epochTime, grand_average_dnone, 'Color', cfg.plotColor{5}, cfg.plotOption{:}); 
xlabel('Time (s)')
ylabel('Amplitude (\muV)');
xline(0);
yline(0);
xlim([-0.1 0.55]); 
ylim([-10 10])
legend('contra low RT', 'ipsi low RT','contra high RT', 'ipsi high RT','none');


ax = gobjects(3,1)
figure('Name', 'grand-average ERP Distractor Left + Right by RT');
ax(1)=subplot(3, 1, 1);
plot(cfg.epochTime, grand_average_contra_lRT, 'Color', cfg.plotColor{1}, cfg.plotOption{:}, 'LineStyle','-'); hold on; 
plot(cfg.epochTime, grand_average_ipsi_lRT, 'Color', cfg.plotColor{2}, cfg.plotOption{:}, 'LineStyle','-'); hold on; 
plot(cfg.epochTime, grand_average_dnone, 'Color', cfg.plotColor{5}, cfg.plotOption{:}); 
legend('contra', 'ipsi','none', 'Location','eastoutside');
title('low RT');
ax(2)=subplot(3, 1, 2);
plot(cfg.epochTime, grand_average_contra_hRT, 'Color', cfg.plotColor{3}, cfg.plotOption{:},'LineStyle','-'); hold on;
plot(cfg.epochTime, grand_average_ipsi_hRT, 'Color', cfg.plotColor{4}, cfg.plotOption{:},'LineStyle','-');
plot(cfg.epochTime, grand_average_dnone, 'Color', cfg.plotColor{5}, cfg.plotOption{:}); 
legend('contra', 'ipsi','none','Location','eastoutside');
title('high RT');
ax(3)=subplot(3, 1, 3);
plot(cfg.epochTime, grand_average_contra_lRT, 'Color', cfg.plotColor{1}, cfg.plotOption{:}, 'LineStyle','-'); hold on; 
plot(cfg.epochTime, grand_average_ipsi_lRT, 'Color', cfg.plotColor{2}, cfg.plotOption{:}, 'LineStyle','-'); hold on; 
plot(cfg.epochTime, grand_average_contra_hRT, 'Color', cfg.plotColor{3}, cfg.plotOption{:},'LineStyle','-'); hold on;
plot(cfg.epochTime, grand_average_ipsi_hRT, 'Color', cfg.plotColor{4}, cfg.plotOption{:},'LineStyle','-');  
plot(cfg.epochTime, grand_average_dnone, 'Color', cfg.plotColor{5}, cfg.plotOption{:}); 
legend('contra low RT', 'ipsi low RT','contra high RT', 'ipsi high RT','no distractor','Location','eastoutside');
title('low + high RT');
xlabel('Time (s)')
for i = 1:length(ax)
    axes(ax(i)); 
    ylabel('Amplitude (\muV)');
    xline(0); 
    yline(0);
    xlim([-0.1 0.55]); 
    ylim([-10 10])
end


%%
%Plot with electrode labels
ax = gobjects(3,1)
figure('Name', 'grand-average ERP');
ax(1)=subplot(3, 1, 1);
plot(cfg.epochTime, PO8_data_dleft, 'Color', cfg.plotColor{1}, cfg.plotOption{:}); hold on;
plot(cfg.epochTime, PO7_data_dleft, 'Color', cfg.plotColor{2}, cfg.plotOption{:}); hold on;
legend('PO8', 'PO7');
title('Distractor on the Left');
ax(2)=subplot(3, 1, 2);
plot(cfg.epochTime, PO7_data_dright, 'Color', cfg.plotColor{1}, cfg.plotOption{:}); hold on;
plot(cfg.epochTime, PO8_data_dright, 'Color', cfg.plotColor{2}, cfg.plotOption{:}); 
legend('PO7', 'PO8');
title('Distractor on the Right');
ax(3)=subplot(3, 1, 3);
plot(cfg.epochTime, grand_average_contra, 'Color', cfg.plotColor{1}, cfg.plotOption{:}); hold on;
plot(cfg.epochTime, grand_average_ipsi, 'Color', cfg.plotColor{2}, cfg.plotOption{:}); 
legend('contra', 'ipsi');
title('Average Left and Right');
xlabel('Time (s)')
for i = 1:length(ax)
    axes(ax(i)); 
    ylabel('Amplitude (\muV)');
    xline(0); 
    yline(0);
    xlim([-0.1 0.55]); 
    ylim([-10 10])
end
%%
% Plot the difference
ax = gobjects(3,1)
figure('Name','Difference Waveform')
ax(1)=subplot(4, 1, 1);
plot(cfg.epochTime, dleft_channel_difference, 'Color', cfg.plotColor{1}, cfg.plotOption{:}); 
title('Distractor on the Left')
ax(2)=subplot(4, 1, 2);
plot(cfg.epochTime, dright_channel_difference, 'Color', cfg.plotColor{1}, cfg.plotOption{:});
title('Distractor on the Right')
ax(3)=subplot(4, 1, 3);
plot(cfg.epochTime, grand_average_channel_difference, 'Color', cfg.plotColor{1}, cfg.plotOption{:});hold on;
xlabel('Time (s)')
for i = 1:length(ax)
    axes(ax(i)); 
    ylabel('Amplitude (\muV)');
    xline(0);
    xline(RT, '--', 'LineWidth', 2); 
    yline(0);
    xlim([-0.1 0.8]); 
    ylim([-20 20])
end
subplot(4, 1, 4);
plot(cfg.epochTime, grand_average_channel_difference, 'Color', cfg.plotColor{1}, cfg.plotOption{:});hold on;
xlabel('Time (s)')
ylabel('Amplitude (\muV)');
xline(0); 
yline(0);
xline(RT, '--', 'LineWidth', 2); 
xlim([-0.1 0.8]); 
ylim([-5 5])
%%
%plot compared to no distractor
cfg.plotColor = {
    [0, 0.4470, 0.7410],      % Blue
    [0.8500, 0.3250, 0.0980], % Orange
    [0.9290, 0.6940, 0.1250], % Yellow
};
figure('Name','ERP Distractor vs No Distractor')
plot(cfg.epochTime, grand_average_contra, 'Color', cfg.plotColor{2}, cfg.plotOption{:}); hold on;
plot(cfg.epochTime, grand_average_ipsi, 'Color', cfg.plotColor{1}, cfg.plotOption{:}); hold on;
plot(cfg.epochTime, grand_average_dnone, 'Color', cfg.plotColor{3}, cfg.plotOption{:}); 
xlabel('Time (s)','FontName', 'Times New Roman')
ylabel('Amplitude (\muV)','FontName', 'Times New Roman');
xline(0);
xline(RT, '--', 'LineWidth', 2); 
yline(0);
xlim([-0.1 0.8]); 
ylim([-6 6])
legend('contra', 'ipsi','none','FontName', 'Times New Roman');
set(gca, 'FontName', 'Times New Roman');

%% topoplot
window_duration = 0.5;  % Time window half-length in seconds
time_points = [0, 0.15, 0.25, 0.35];  % Time points of interest
srate=512
fixed_color_limits = [-5, 5];
channels_to_plot = {cfg.chanlocs.labels}
plot_channel_indices = find(ismember({cfg.chanlocs.labels}, channels_to_plot));

for i = 1:length(time_points)
    current_time = time_points(i);
    window_start_time = current_time - window_duration;
    window_end_time = current_time + window_duration;

    window_start_time = max(window_start_time, cfg.epochTime(1));
    window_end_time = min(window_end_time, cfg.epochTime(end));

    % Find the closest indices in cfg.epochTime
    [~, window_start_idx] = min(abs(cfg.epochTime - window_start_time));
    [~, window_end_idx] = min(abs(cfg.epochTime - window_end_time));

    %avg amp across time
    avg_distractor = mean(epochs.data(window_start_idx:window_end_idx, :, epochs.labels == 2), 1);
    avg_no_distractor = mean(epochs.data(window_start_idx:window_end_idx, :, epochs.labels == 0), 1);

    % average amplitude across trials
    avg_distractor_amp = mean(avg_distractor, 3);
    avg_no_distractor_amp = mean(avg_no_distractor, 3);


    plot_chanlocs = cfg.chanlocs(plot_channel_indices);

    figure;
    subplot(1, 2, 1);
    topoplot(avg_distractor_amp, plot_chanlocs,'maplimits', fixed_color_limits,'electrodes', 'labels');
    colorbar;
    title(sprintf('Distractor, Time = %.2f s', time_points(i)));

    subplot(1, 2, 2);
    topoplot(avg_no_distractor_amp,plot_chanlocs,'maplimits', fixed_color_limits,'electrodes', 'labels');
    colorbar;
    title(sprintf('No Distractor, Time = %.2f s', time_points(i)));
end



%%
% channelsToInclude = setdiff(1:size(epochs.data, 2), 32)
% plot(cfg.epochTime, mean(epochs.data(:, channelsToInclude, epochs.labels == true), 3), 'Color', cfg.plotColor{1}, cfg.plotOption{:}); hold on;
% 
% 
% figure('Name', 'all-channels');
% gap = 2;
% plot(cfg.epochTime, bsxfun(@plus, mean(epochs.data(:, :, epochs.labels == false), 3), 0:gap:gap*(size(epochs.data,2)-1)), 'Color', cfg.plotColor{1}, cfg.plotOption{:}); hold on;
% plot(cfg.epochTime, bsxfun(@plus, mean(epochs.data(:, :, epochs.labels == true), 3), 0:gap:gap*(size(epochs.data,2)-1)), 'Color', cfg.plotColor{2}, cfg.plotOption{:});

figure('Name', 'single-trial distractor either left or right')
subplot(3, 1, 1);
imagesc(cfg.epochTime, 1:sum(epochs.labels == 0), squeeze(epochs.data(:, 6, epochs.labels==0))'); caxis([-20 20]);
title('No Distractor')
xlim([-0.1 0.4])
subplot(3, 1, 2);
imagesc(cfg.epochTime, 1:sum((epochs.labels == 2)), squeeze(epochs.data(:, 6, epochs.labels==2))');  caxis([-20 20]);
title('Distractor on the Left')
xlim([-0.1 0.4])
subplot(3, 1, 3);
imagesc(cfg.epochTime, 1:sum(epochs.labels == 1), squeeze(epochs.data(:, 6, epochs.labels==1))'); caxis([-20 20]);
title('Distractor on the Right')
xlim([-0.1 0.4])



figure('Name', 'single-trial distractor either left or right')
subplot(4, 1, 1);
imagesc(cfg.epochTime, 1:sum(epochs.labels == 0), squeeze(epochs.data(:, PO7, epochs.labels==0))'); caxis([-20 20]);
title('No Distractor PO7')
xlim([-0.1 0.7])
subplot(4, 1, 2);
imagesc(cfg.epochTime, 1:sum((epochs.labels == 2)), squeeze(epochs.data(:, PO7, epochs.labels==1))');  caxis([-20 20]);
title('Distractor on the Left PO7')
xlim([-0.1 0.7])
subplot(4, 1, 3);
imagesc(cfg.epochTime, 1:sum(epochs.labels == 0), squeeze(epochs.data(:, PO8, epochs.labels==0))'); caxis([-20 20]);
title('No Distractor PO8')
xlim([-0.1 0.7])
subplot(4, 1, 4);
imagesc(cfg.epochTime, 1:sum((epochs.labels == 2)), squeeze(epochs.data(:, PO8, epochs.labels==1))');  caxis([-20 20]);
title('Distractor on the Left PO8')
xlim([-0.1 0.7])

figure('Name', 'single-trial')
subplot(4, 1, 1);
imagesc(cfg.epochTime, 1:sum(((epochs.labels == 0))), squeeze(epochs.data(:, PO7, ((epochs.labels == 0))))'); caxis([-20 20]);
title('No Distractor PO7')
%xlim([0.1 0.3])
subplot(4, 1, 2);
imagesc(cfg.epochTime, 1:sum(((epochs.labels == 1))), squeeze(epochs.data(:, PO7, ((epochs.labels == 1))))');  caxis([-20 20]);
title('Distractor on the Right PO7')
%xlim([0.1 0.3])
subplot(4, 1, 3);
imagesc(cfg.epochTime, 1:sum(((epochs.labels == 0))), squeeze(epochs.data(:, PO8, ((epochs.labels == 0))))'); caxis([-20 20]);
title('No Distractor PO8')
%xlim([0.1 0.3])
subplot(4, 1, 4);
imagesc(cfg.epochTime, 1:sum(((epochs.labels == 1))), squeeze(epochs.data(:, PO8, ((epochs.labels == 1))))');  caxis([-20 20]);
title('Distractor on the Right PO8')
%xlim([0.1 0.3])




% figure('Name', 'single-trial')
% subplot(2, 1, 1);
% imagesc(cfg.epochTime, 1:sum(epochs.labels == 0), squeeze(singletrial_difference_dleft)'); caxis([-20 20]);
% title('No Distractor contra-ipsi')
% xlim([-0.1 0.7])
% subplot(2, 1, 2);
% imagesc(cfg.epochTime, 1:sum(epochs.labels == 1), squeeze(singletrial_difference_dleft)');  caxis([-20 20]);
% title('Distractor on the Left contra-ipsi')
% xlim([-0.1 0.7])


figure('Name', 'single-trial')
subplot(2, 1, 1);
imagesc(cfg.epochTime, 1:sum(epochs.labels == 0), squeeze(singletrial_difference_dright)'); caxis([-20 20]);
title('No Distractor contra-ipsi')
xlim([-0.1 0.7])
subplot(2, 1, 2);
imagesc(cfg.epochTime, 1:sum(epochs.labels == 1), squeeze(singletrial_difference_dright)');  caxis([-20 20]);
title('Distractor on the Right contra-ipsi')
xlim([-0.1 0.7])


figure('Name', 'r-squared')
r2 = compute_r2(epochs.data, epochs.labels);
imagesc(cfg.epochTime, 1:size(epochs.data, 2), r2');
xlim([-0.1 0.7])

drawnow;

%%%%%%%%%%%%%%%%%%%%
%% Classification %%
%%%%%%%%%%%%%%%%%%%%
%quadratic diagLinear diagQuadratic pseudoLinear pseudoQuadratic
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
end