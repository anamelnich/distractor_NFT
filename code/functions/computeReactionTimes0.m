function training = computeReactionTimes(training, index, sampling_rate)
    % Computes reaction times and labels based on trigger information.
    %
    % Inputs:
    %   training (struct): Training data structure containing 'trigger' field.
    %   set_size (int): Set size of the trial (e.g., 4, 6, 10).
    %   index (struct): Index structure with 'pos' and 'typ' fields.
    %   sampling_rate (double): Sampling rate in Hz (e.g., 512).
    %
    % Outputs:
    %   training (struct): Updated training structure with RT.time, RT.condition, and RT.typ fields.
    
    % Validate inputs
    if ~isfield(training, 'trigger')
        error('The training structure must contain a ''trigger'' field.');
    end
    
    if nargin < 4 || isempty(sampling_rate)
        sampling_rate = 512;  % Default sampling rate
        warning('Sampling rate not provided. Using default value of 512 Hz.');
    end
    
    triggers = training.trigger;  % Trigger codes

    trial_start_positions = index.pos;  % Positions of trial start triggers
    trial_types = index.typ;  % Types of trial start triggers (1, 2, 0)
    
    % Define response trigger codes
    response_triggers = [11, 12, 13, 21, 22, 23];
    
    % Initialize arrays to store reaction times and conditions
    num_trials = length(trial_start_positions);
    reaction_times = NaN(num_trials, 1);  % Preallocate with NaN
    conditions = zeros(num_trials, 1);  % To store condition codes
    
    % Loop through each trial to find response and calculate RT
    for i = 1:num_trials
        start_pos = trial_start_positions(i);
        
        % Ensure we don't exceed array bounds
        if start_pos >= length(triggers)
            warning('Trial start position %d exceeds trigger array length.', start_pos);
            continue;
        end
        
        % Find the next response trigger after the trial start
        response_relative_pos = find(ismember(triggers(start_pos+1:end), response_triggers), 1, 'first');
        
        if ~isempty(response_relative_pos)
            response_pos = start_pos + response_relative_pos;  % Absolute position in trigger array
            % Calculate reaction time in seconds
            rt = (response_pos - start_pos) / sampling_rate;
            reaction_times(i) = rt;
            
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
    
    % Store RTs and Conditions in the training structure
    training.RT.time = reaction_times;
    training.RT.condition = conditions;
    
    % Compute average RT
    average_rt = mean(reaction_times);
    
    % Assign labels: 1 for below average RT, 0 for above or equal
    labels = double(reaction_times < average_rt);
    training.RT.typ = labels;
    
    % Display Statistics
    fprintf('Computed Reaction Times for %d trials.\n', length(reaction_times));
    fprintf('Average Reaction Time: %.3f seconds\n', average_rt);
    fprintf('Trials Below Average RT: %d\n', sum(labels));
    fprintf('Trials Above or Equal to Average RT: %d\n', sum(~labels));
    
    % Optional: Display a summary table
    % training.RT.table = table(trial_start_positions, reaction_times, labels, ...
    %                           'VariableNames', {'TrialStartPos', 'ReactionTime', 'Label'});
end






