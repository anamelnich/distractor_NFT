function training = computeReactionTimes(training, index, sampling_rate)
    % Computes reaction times and labels based on trigger information.
    %
    % Inputs:
    %   training (struct): Training data structure containing 'trigger' field.
    %   index (struct): Index structure with 'pos' and 'typ' fields.
    %   sampling_rate (double): Sampling rate in Hz (e.g., 512).
    %
    % Outputs:
    %   training (struct): Updated training structure with RT.time, RT.condition, RT.typ, and RT.resp fields.

    % Validate inputs
    if ~isfield(training, 'trigger')
        error('The training structure must contain a ''trigger'' field.');
    end

    if nargin < 3 || isempty(sampling_rate)
        sampling_rate = 512;  % Default sampling rate
        warning('Sampling rate not provided. Using default value of 512 Hz.');
    end

    triggers = training.trigger;  % Trigger codes

    trial_start_positions = index.pos;  % Positions of trial start triggers
    trial_types = index.typ;  % Types of trial start triggers (e.g., 1, 2, 0)

    % Define response trigger codes
    response_triggers = [11, 12, 13, 21, 22, 23];

    % Initialize arrays to store reaction times, conditions, and response codes
    num_trials = length(trial_start_positions);
    reaction_times = NaN(num_trials, 1);  % Preallocate with NaN
    conditions = zeros(num_trials, 1);  % To store condition codes
    response_codes = NaN(num_trials, 1);  % To store response trigger codes

    % Loop through each trial to find response and calculate RT
    for i = 1:num_trials
        start_pos = trial_start_positions(i);

        % Determine the end position for the current trial
        if i < num_trials
            end_pos = trial_start_positions(i + 1) - 1;  % Up to the next trial start
        else
            end_pos = length(triggers);  % Last trial goes to the end of the triggers
        end

        % Ensure we don't exceed array bounds
        if start_pos >= length(triggers)
            warning('Trial start position %d exceeds trigger array length.', start_pos);
            continue;
        end

        % Search for response triggers between start_pos and end_pos
        triggers_in_trial = triggers(start_pos+1:end_pos);

        % Find the first occurrence of a response trigger
        response_idx = find(ismember(triggers_in_trial, response_triggers), 1, 'first');

        if ~isempty(response_idx)
            response_pos = start_pos + response_idx;
            % Calculate reaction time in seconds
            rt = (response_pos - start_pos) / sampling_rate;
            reaction_times(i) = rt;

            % Store the response trigger code
            response_code = triggers(response_pos);
            response_codes(i) = response_code;

            % Assign condition based on trial type
            conditions(i) = trial_types(i);
        else
            % No response found; RT remains NaN
            warning('No response trigger found for trial starting at position %d.', start_pos);
        end
    end

    % Remove trials with NaN RTs (no response found)
    valid_trials = ~isnan(reaction_times);
    reaction_times = reaction_times(valid_trials);
    conditions = conditions(valid_trials);
    trial_start_positions = trial_start_positions(valid_trials);
    response_codes = response_codes(valid_trials);

    % Store RTs and Conditions in the training structure
    training.RT.time = reaction_times;
    training.RT.condition = conditions;

    % Compute average RT
    average_rt = mean(reaction_times);

    % Assign labels for RT.typ based on reaction time
    labels = zeros(length(reaction_times), 1);  % Initialize labels

    labels(reaction_times < average_rt) = 1;
    labels(reaction_times >= average_rt) = 0;

    % Store labels in training structure
    training.RT.typ = labels;

    % Assign RT.resp based on response codes
    resp_labels = NaN(length(response_codes), 1);  % Initialize resp_labels

    % Assign resp_labels: 1 if response code is 11 or 21, 0 if response code is 12,13,22,23
    resp_labels(ismember(response_codes, [11, 21])) = 1;
    resp_labels(ismember(response_codes, [12, 13, 22, 23])) = 0;

    % Store resp_labels in training structure
    training.RT.resp = resp_labels;

    % Display Statistics
    fprintf('Computed Reaction Times for %d trials.\n', length(reaction_times));
    fprintf('Average Reaction Time: %.3f seconds\n', average_rt);
    fprintf('Trials with RT < average (RT.typ=1): %d\n', sum(labels == 1));
    fprintf('Trials with RT >= average (RT.typ=0): %d\n', sum(labels == 0));
    fprintf('Trials with response code 11 or 21 (RT.resp=1): %d\n', sum(resp_labels == 1));
    fprintf('Trials with response code 12, 13, 22, or 23 (RT.resp=0): %d\n', sum(resp_labels == 0));

    % Optional: Display a summary table
    % training.RT.table = table(trial_start_positions, reaction_times, labels, resp_labels, response_codes, ...
    %                           'VariableNames', {'TrialStartPos', 'ReactionTime', 'RT_typ', 'RT_resp', 'ResponseCode'});
end




