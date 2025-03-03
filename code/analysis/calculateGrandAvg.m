function averages = calculateGrandAvg(data, labels, leftElectrodeIndices, rightElectrodeIndices)
%calculateGrandAvg(trainEpochs, trainLabels, leftElectrodeIndices, rightElectrodeIndices);
% calculateGrandAvg computes trial-level averages, grand averages, and standard
% deviations for each condition. It also computes a difference wave for the
% no-distractor condition with randomized subtraction.
%
% Inputs:
%   data - 3D array [time_points x electrodes x trials]
%   labels - vector indicating condition for each trial
%   leftElectrodeIndices - indices for left electrodes
%   rightElectrodeIndices - indices for right electrodes
%
% Output:
%   averages - structure containing trial-level data, grand averages, and SDs

%% Condition 1: Distractor Right (label == 1)
idx_dright = (labels == 1);
% Average over electrodes to get trial-level data [time x nTrials_dright]
trial_leftElect_dright = squeeze(mean(data(:, leftElectrodeIndices, idx_dright), 2));
trial_rightElect_dright = squeeze(mean(data(:, rightElectrodeIndices, idx_dright), 2));

% Compute grand average (across trials) and standard deviation (per time point)
grand_avg_leftElect_dright = mean(trial_leftElect_dright, 2);
std_leftElect_dright     = std(trial_leftElect_dright, 0, 2);
grand_avg_rightElect_dright = mean(trial_rightElect_dright, 2);
std_rightElect_dright     = std(trial_rightElect_dright, 0, 2);

%% Condition 2: Distractor Left (label == 2)
idx_dleft = (labels == 2);
trial_rightElect_dleft = squeeze(mean(data(:, rightElectrodeIndices, idx_dleft), 2));
trial_leftElect_dleft  = squeeze(mean(data(:, leftElectrodeIndices, idx_dleft), 2));

grand_avg_rightElect_dleft = mean(trial_rightElect_dleft, 2);
std_rightElect_dleft     = std(trial_rightElect_dleft, 0, 2);
grand_avg_leftElect_dleft  = mean(trial_leftElect_dleft, 2);
std_leftElect_dleft      = std(trial_leftElect_dleft, 0, 2);

%% Condition 0: No Distractor (label == 0)
idx_dnone = (labels == 0);
trial_leftElect_dnone  = squeeze(mean(data(:, leftElectrodeIndices, idx_dnone), 2));
trial_rightElect_dnone = squeeze(mean(data(:, rightElectrodeIndices, idx_dnone), 2));

grand_avg_leftElect_dnone  = mean(trial_leftElect_dnone, 2);
std_leftElect_dnone      = std(trial_leftElect_dnone, 0, 2);
grand_avg_rightElect_dnone = mean(trial_rightElect_dnone, 2);
std_rightElect_dnone     = std(trial_rightElect_dnone, 0, 2);

%% Combined Grand Averages for Distractor Conditions
% For distractor left: use right electrodes as contralateral and left electrodes as ipsilateral.
grand_average_contra_dleft = grand_avg_rightElect_dleft;
grand_average_ipsi_dleft   = grand_avg_leftElect_dleft;
% For distractor right: use left electrodes as contralateral and right electrodes as ipsilateral.
grand_average_contra_dright = grand_avg_leftElect_dright;
grand_average_ipsi_dright   = grand_avg_rightElect_dright;

% Combined contra and ipsi averages
combined_trials_contra = [trial_rightElect_dleft, trial_leftElect_dright];
combined_trials_ipsi = [trial_leftElect_dleft, trial_rightElect_dright];

grand_average_contra = mean(combined_trials_contra,2);
std_grand_average_contra = std(combined_trials_contra,0,2);

grand_average_ipsi = mean(combined_trials_ipsi,2);
std_grand_average_ipsi = std(combined_trials_ipsi,0,2);
 
%% Combined Grand Average for No Distractor
combined_trials_dnone = [trial_rightElect_dnone, trial_leftElect_dnone];

grand_average_dnone = mean(combined_trials_dnone,2);
std_grand_average_dnone = std(combined_trials_dnone,0,2);

%% Distractor Difference Wave (Randomized Subtraction)

trial_d_diff = combined_trials_contra - combined_trials_ipsi;
grand_average_d_diff = mean(trial_d_diff, 2);
std_d_diff = std(trial_d_diff, 0, 2);

%% No Distractor Difference Wave (Randomized Subtraction)
nTrials_dnone = sum(idx_dnone);
% Generate a random sign (+1 or -1) for each trial.
signs = (rand(1, nTrials_dnone) > 0.5) * 2 - 1;
% For each trial, if sign==+1, compute (left - right); if sign==-1, compute (right - left).
trial_dnone_diff = signs .* (trial_leftElect_dnone - trial_rightElect_dnone);
grand_average_dnone_diff = mean(trial_dnone_diff, 2);
std_dnone_diff = std(trial_dnone_diff, 0, 2);

%% Store results

averages.grand_average_contra_dleft      = grand_average_contra_dleft;
averages.grand_average_contra_dleft_sd   = std_rightElect_dleft;
averages.grand_average_contra_dright     = grand_average_contra_dright;
averages.grand_average_contra_dright_sd  = std_leftElect_dright;
averages.grand_average_ipsi_dleft        = grand_average_ipsi_dleft;
averages.grand_average_ipsi_dleft_sd     = std_leftElect_dleft;
averages.grand_average_ipsi_dright       = grand_average_ipsi_dright;
averages.grand_average_ipsi_dright_sd    = std_rightElect_dright;

averages.grand_average_contra = grand_average_contra;
averages.grand_average_contra_sd = std_grand_average_contra;
averages.grand_average_ipsi   = grand_average_ipsi;
averages.grand_average_ipsi_sd = std_grand_average_ipsi;

averages.grand_average_dnone_left   = grand_avg_leftElect_dnone;
averages.grand_average_dnone_left_sd = std_leftElect_dnone;
averages.grand_average_dnone_right   = grand_avg_rightElect_dnone;
averages.grand_average_dnone_right_sd = std_rightElect_dnone;

averages.grand_average_dnone = grand_average_dnone;
averages.grand_average_dnone_sd = std_grand_average_dnone;

% Distractor Difference Wave:
averages.grand_average_d_diff      = grand_average_d_diff;
averages.grand_average_d_diff_sd   = std_d_diff;
averages.trial_d_diff              = trial_d_diff;

% No Distractor Difference Wave:
averages.grand_average_dnone_diff      = grand_average_dnone_diff;
averages.grand_average_dnone_diff_sd   = std_dnone_diff;
averages.trial_dnone_diff              = trial_dnone_diff;

end
