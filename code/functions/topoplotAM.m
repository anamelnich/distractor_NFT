% --- Configuration Section ---
window_duration = 0.05;  % Time window half-length in seconds (total window = 0.1 seconds)
time_points = [0, 0.15, 0.25, 0.3];  % Time points of interest in seconds
srate = 512;  % Sampling rate in Hz
fixed_color_limits = [-2, 2];  % Color limits for topoplot

% Define channels to plot
% Specify the channels you want to include. Example: {'PO7', 'PO8'}
channels_to_plot = {cfg.chanlocs.labels};  % Adjust as needed

% Find channel indices
plot_channel_indices = find(ismember({cfg.chanlocs.labels}, channels_to_plot));

% Validate channel selection
if isempty(plot_channel_indices)
    error('None of the specified channels found in cfg.chanlocs.');
end

% Retrieve channel locations once
plot_chanlocs = cfg.chanlocs(plot_channel_indices);

% Validate cfg.epochTime
if ~isfield(cfg, 'epochTime')
    error('cfg.epochTime is not defined.');
end

% --- Plotting Loop ---
for i = 1:length(time_points)
    current_time = time_points(i);
    
    % Define the start and end times of the window in seconds
    window_start_time = current_time - window_duration;
    window_end_time = current_time + window_duration;
    
    % Ensure the window doesn't exceed epoch boundaries
    window_start_time = max(window_start_time, cfg.epochTime(1));
    window_end_time = min(window_end_time, cfg.epochTime(end));
    
    % Find the closest indices in cfg.epochTime for the window
    [~, window_start_idx] = min(abs(cfg.epochTime - window_start_time));
    [~, window_end_idx] = min(abs(cfg.epochTime - window_end_time));
    
    % Validate window indices
    if window_start_idx > window_end_idx
        warning('For time point %.2f s, window_start_idx (%d) > window_end_idx (%d). Skipping this window.', ...
            current_time, window_start_idx, window_end_idx);
        continue;  % Skip to the next iteration
    end
    
    % Define trial conditions
    distractor_right_trials = (epochs.labels == 1) %& (epochs.RT == 0);
    %distractor_left_trials = (epochs.labels == 2) %& (epochs.RT == 0);
    no_distractor_trials = (epochs.labels == 0) %& (epochs.RT == 0);
    
    % Check if there are trials for each condition
    if ~any(distractor_right_trials)
        warning('No Distractor Right trials found at time point %.2f s. Skipping Distractor Right plot.', current_time);
    end
    % if ~any(distractor_left_trials)
    %     warning('No Distractor Left trials found at time point %.2f s. Skipping Distractor Left plot.', current_time);
    % end
    if ~any(no_distractor_trials)
        warning('No No Distractor trials found at time point %.2f s. Skipping No Distractor plot.', current_time);
    end
    
    % Extract and average data for Distractor Right (Label 2)
    if any(distractor_right_trials)
        data_distractor_right = epochs.data(window_start_idx:window_end_idx, plot_channel_indices, distractor_right_trials);
        avg_distractor_right = squeeze(mean(mean(data_distractor_right, 1), 3));  % [channels x 1]
    else
        avg_distractor_right = NaN(length(plot_channel_indices), 1);
    end
    
    % Extract and average data for Distractor Left (Label 1)
    % if any(distractor_left_trials)
    %     data_distractor_left = epochs.data(window_start_idx:window_end_idx, plot_channel_indices, distractor_left_trials);
    %     avg_distractor_left = squeeze(mean(mean(data_distractor_left, 1), 3));  % [channels x 1]
    % else
    %     avg_distractor_left = NaN(length(plot_channel_indices), 1);
    % end
    
    % Extract and average data for No Distractor (Label 0)
    if any(no_distractor_trials)
        data_no_distractor = epochs.data(window_start_idx:window_end_idx, plot_channel_indices, no_distractor_trials);
        avg_no_distractor = squeeze(mean(mean(data_no_distractor, 1), 3));  % [channels x 1]
    else
        avg_no_distractor = NaN(length(plot_channel_indices), 1);
    end
    
    % --- Plotting Section ---
    figure('Name', sprintf('Grand-average ERP at %.2f s', current_time));
    
    % Create a 1x3 subplot layout for the three conditions
    subplot(1, 2, 1);
    if ~all(isnan(avg_distractor_right))
        topoplot(avg_distractor_right, plot_chanlocs, 'maplimits', fixed_color_limits, ...
            'electrodes', 'labels');
        colorbar;
        title(sprintf('Distractor \nTime = %.2f s', current_time), 'FontSize', 10);
    else
        text(0.5, 0.5, 'No Data', 'HorizontalAlignment', 'center', 'FontSize', 12);
        title(sprintf('Distractor Right\nTime = %.2f s', current_time), 'FontSize', 10);
    end
    
    % subplot(1, 3, 2);
    % if ~all(isnan(avg_distractor_left))
    %     topoplot(avg_distractor_left, plot_chanlocs, 'maplimits', fixed_color_limits, ...
    %         'electrodes', 'labels');
    %     colorbar;
    %     title(sprintf('Distractor Left\nTime = %.2f s', current_time), 'FontSize', 10);
    % else
    %     text(0.5, 0.5, 'No Data', 'HorizontalAlignment', 'center', 'FontSize', 12);
    %     title(sprintf('Distractor Left\nTime = %.2f s', current_time), 'FontSize', 10);
    % end
    
    subplot(1, 2, 2);
    if ~all(isnan(avg_no_distractor))
        topoplot(avg_no_distractor, plot_chanlocs, 'maplimits', fixed_color_limits, ...
            'electrodes', 'labels');
        colorbar;
        title(sprintf('No Distractor\nTime = %.2f s', current_time), 'FontSize', 10);
    else
        text(0.5, 0.5, 'No Data', 'HorizontalAlignment', 'center', 'FontSize', 12);
        title(sprintf('No Distractor\nTime = %.2f s', current_time), 'FontSize', 10);
    end
    
    % Overall Figure Title (requires MATLAB R2018b or later)
    if verLessThan('matlab', '9.4')  % R2018b is version 9.4
        suptitle(sprintf('Grand-average ERP at Time = %.2f s', current_time));
    else
        sgtitle(sprintf('Grand-average ERP at Time = %.2f s', current_time));
    end
end

