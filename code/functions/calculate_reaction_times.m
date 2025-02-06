function rt_table = calculate_reaction_times(trigger_file)
    % Calculates reaction times from combined EEG data.
    %
    % Inputs:
    %   trigger_file (matrix): Combined data with columns [Trial, TriggerCode, Time].
    %
    % Outputs:
    %   rt_table (table): Table containing Trial, Condition, and ReactionTime.
    
    % Extract columns
    trials = trigger_file(:, 1);
    trigger_codes = trigger_file(:, 2);
    times = trigger_file(:, 3);
    
    % Initialize cell array for trigger codes as strings
    trigger_codes_str = arrayfun(@num2str, trigger_codes, 'UniformOutput', false);
    
    % Initialize variables to store reaction times
    reaction_times = [];
    
    % Get unique trial numbers
    unique_trials = unique(trials);
    
    for i = 1:length(unique_trials)
        trial = unique_trials(i);
        % Indices for the current trial
        trial_indices = find(trials == trial);
        
        % Get trigger codes and times for the current trial
        trial_triggers = trigger_codes_str(trial_indices);
        trial_times = times(trial_indices);
        
        % Identify trial start trigger (3-digit starting with '1' or '2')
        trial_start_indices = find(cellfun(@(x) length(x) == 3 && (x(1) == '1' || x(1) == '2'), trial_triggers));
        
        % Identify response trigger (2-digit starting with '1' or '2')
        response_indices = find(cellfun(@(x) length(x) == 2 && (x(1) == '1' || x(1) == '2'), trial_triggers));
        
        if ~isempty(trial_start_indices) && ~isempty(response_indices)
            % Assume the first occurrence is the trial start and response
            trial_start_idx = trial_start_indices(1);
            response_idx = response_indices(1);
            
            % Ensure response occurs after trial start
            if response_idx > trial_start_idx
                rt = trial_times(response_idx) - trial_times(trial_start_idx);
                condition = str2double(trial_triggers{trial_start_idx});
                
                % Append to reaction times list
                reaction_times(end+1, :) = [trial, condition, rt];
            end
        end
    end
    
    % Convert to table for better readability
    rt_table = array2table(reaction_times, 'VariableNames', {'Trial', 'Condition', 'ReactionTime'});
end
