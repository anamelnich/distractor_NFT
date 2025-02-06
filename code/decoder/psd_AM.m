% Assuming epochs.data is [time x channels x trials]
% Extract the indices for each condition
idx_no_distractor = epochs.labels == 0;
idx_distractor = epochs.labels ~= 0;  % labels 1 and 2 represent distractors

% Extract data for each condition
data_no_distractor = epochs.data(:, :, idx_no_distractor);
data_distractor = epochs.data(:, :, idx_distractor);

% Sampling frequency
fs = cfg.fsamp;

% Time indices for PSD computation (ensure this is correctly defined)
timeIndices = cfg.psd.time;

% PSD parameters
window = cfg.psd.window;
nfft = cfg.psd.nfft;
noverlap = cfg.psd.overlap;

% If overlap is not specified, set a default value (e.g., 50% overlap)
if isempty(noverlap)
    noverlap = floor(length(window) / 2);
end

% Initialize variables to store PSDs
psds_no = [];
psds_distractor = [];

% Compute PSDs for No Distractor condition
for trial = 1:size(data_no_distractor, 3)
    data_trial = mean(data_no_distractor(timeIndices, :, trial), 2);  % Average over channels
    [Pxx, f] = pwelch(data_trial, window, noverlap, nfft, fs);
    psds_no(:, trial) = Pxx;
end
psd_no_mean = mean(psds_no, 2);

% Compute PSDs for Distractor condition
for trial = 1:size(data_distractor, 3)
    data_trial = mean(data_distractor(timeIndices, :, trial), 2);  % Average over channels
    [Pxx, ~] = pwelch(data_trial, window, noverlap, nfft, fs);
    psds_distractor(:, trial) = Pxx;
end
psd_distractor_mean = mean(psds_distractor, 2);

% Plot the PSDs
figure;
plot(f, 10*log10(psd_no_mean), 'b', 'LineWidth', 2);
hold on;
plot(f, 10*log10(psd_distractor_mean), 'r', 'LineWidth', 2);
xlabel('Frequency (Hz)');
ylabel('Power/Frequency (dB/Hz)');
legend('No Distractor', 'Distractor');
title('Power Spectral Density');
grid on;


