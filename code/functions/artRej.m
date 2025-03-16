function badEpochs = artRej(trainEpochs, params)
% Artifact Rejection Script
% =====================================
% This script performs artifact rejection on EEG epochs based on the following criteria:
% 1. Voltage step per sampling point > 30 μV
% 2. Absolute amplitude exceeding ±60 μV
%
% Inputs:
% - trainEpochs: 3D matrix of EEG data [time x channel x epoch] (768 x 58 x 240)
% - params: Struct containing various parameters, including sampling rate (fsamp)
%
% Outputs:
% - badEpochs: Logical array of bad epochs

% Define Artifact Thresholds
voltageStepThreshold = 50;      
absoluteAmplitudeThreshold = 75; 

timeIdx = params.epochRejection.time;
epochData = trainEpochs(timeIdx, :, :);
% Define Window for Activity Change (500 ms)
windowDurationSec = 0.5;  % 500 ms
windowSamples = round(windowDurationSec * params.fsamp);  

% Dimensions of trainEpochs
[n_time, n_channels, n_epochs] = size(epochData);

% Initialize Logical Arrays to Mark Bad Epochs
badEpoch_VoltageStep = false(1, n_epochs);
badEpoch_AbsAmplitude = false(1, n_epochs);

% =====================================
% 1. Voltage Step Artifact Detection
% =====================================
% 
% % Compute the difference between consecutive samples for all epochs and channels
% % Resulting size: [n_time-1 x n_channels x n_epochs]
% voltageSteps = abs(diff(epochData, 1, 1));
% 
% % Check if any voltage step exceeds the threshold in any channel for each epoch
% % Resulting size after any(): [1 x 1 x n_epochs]
% voltageStepExceeds = any(any(voltageSteps > voltageStepThreshold, 1), 2);  
% 
% % Convert to logical array [1 x n_epochs]
% badEpoch_VoltageStep = squeeze(voltageStepExceeds);


% =====================================
% 2. Absolute Amplitude Artifact Detection
% =====================================

% Check if any sample in any channel exceeds the absolute amplitude threshold
% Resulting size: [n_time x n_channels x n_epochs]
absAmplitudes = abs(epochData);

% Check if any absolute amplitude exceeds the threshold
absAmplitudeExceeds = any(any(absAmplitudes > absoluteAmplitudeThreshold, 1), 2);

% Convert to logical array [1 x n_epochs]
badEpoch_AbsAmplitude = squeeze(absAmplitudeExceeds);

% =====================================
% 3. Combine All Artifact Conditions
% =====================================

% An epoch is bad if any of the three conditions are true
% badEpochs = badEpoch_VoltageStep | badEpoch_AbsAmplitude;
badEpochs = badEpoch_AbsAmplitude;

end