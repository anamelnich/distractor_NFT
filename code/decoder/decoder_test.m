%% Define the electrode names for left and right
LeftElectrodes = {'P1', 'P3', 'P5', 'P7', 'PO3', 'PO7'};
RightElectrodes = {'P2', 'P4', 'P6', 'P8', 'PO4', 'PO8'};

%% Map electrode names to indices
leftElectrodeIndices = find(ismember({cfg.chanlocs.labels},LeftElectrodes));
rightElectrodeIndices = find(ismember({cfg.chanlocs.labels},RightElectrodes));


%% Initialize epochs.data with the correct dimensions
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

    % If label == 2, you might want to rearrange the data to match the order of left electrodes
    if label == 2
        % Optionally, you can rearrange the right electrodes to match the left electrodes
        % This step ensures that each channel corresponds to the same scalp location across trials
        % For example:
        % [P2, P4, P6, P8, PO4, PO8] corresponds to [P1, P3, P5, P7, PO3, PO7]
        % You can create a mapping if needed
    end

    % Store the data in epochs.data
    epochs.data(:, :, i_trial) = epochData;

    % Find the file ID for the current trial
    temp = find(training.index.pos(i_trial) <= training.eof, 1, 'first');
    epochs.file_id(i_trial) = temp;
end

%% decoder

n_files = length(training.eof);
epochs.posteriors = nan(n_trials, 1);
for i_file = 1:n_files
    train_index = epochs.file_id ~= i_file; 
    test_index = epochs.file_id == i_file;
    decoder = computeDecoder(epochs.data(:, :, train_index), epochs.labels(train_index), cfg);
    epochs.posteriors(test_index) = singleClassification(decoder, epochs.data(:, :, test_index));
end

% Continue with performance evaluation and saving the decoder
