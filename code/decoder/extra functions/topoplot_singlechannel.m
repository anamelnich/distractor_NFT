% Define timepoints of interest
window_duration = 0.5;  % Time window half-length in seconds
time_points = [0, 0.15, 0.25, 0.35];  % Time points of interest
srate = 512;

% Define the channels you want to plot
channels_to_plot = {'PO7', 'PO8'}; % Specify your desired channels
plot_channel_indices = find(ismember({cfg.chanlocs.labels}, channels_to_plot));

% Check if all specified channels are found
if length(plot_channel_indices) ~= length(channels_to_plot)
    error('One or more specified channels not found in cfg.chanlocs.');
end

% Define fixed color scale limits
fixed_color_limits = [-2, 2];

% Determine the number of time points and conditions
num_timepoints = length(time_points);
num_conditions = 2; % Distractor and No Distractor

% Create a single figure for all topoplots
figure('Name', 'Topoplots for PO7 and PO8 Across Time Points and Conditions', ...
       'NumberTitle', 'off', ...
       'Position', [100, 100, 1600, 800]);

% Loop through each time point
for i = 1:num_timepoints
    % Define window times in seconds
    window_start_time = time_points(i) - window_duration;
    window_end_time = time_points(i) + window_duration;
    
    % Clamp window times to epoch bounds
    window_start_time = max(window_start_time, cfg.epochTime(1));
    window_end_time = min(window_end_time, cfg.epochTime(end));
    
    % Find the closest sample indices in cfg.epochTime
    [~, window_start_idx] = min(abs(cfg.epochTime - window_start_time));
    [~, window_end_idx] = min(abs(cfg.epochTime - window_end_time));
    
    % Check for valid indices
    if window_start_idx < 1 || window_end_idx > size(epochs.data, 1) || window_start_idx >= window_end_idx
        warning('Invalid window indices for time point %.2f s. Skipping.', time_points(i));
        continue;
    end
    
    % Average amplitude across time within the window for each condition
    avg_distractor = mean(epochs.data(window_start_idx:window_end_idx, :, epochs.labels == 1), 1); % [1 x channels x trials]
    avg_no_distractor = mean(epochs.data(window_start_idx:window_end_idx, :, epochs.labels == 2), 1); % [1 x channels x trials]
    
    % Average amplitude across trials
    avg_distractor_amp = mean(avg_distractor, 3); % [1 x channels]
    avg_no_distractor_amp = mean(avg_no_distractor, 3); % [1 x channels]
    
    % Subset the data and channel locations for plotting
    plot_distractor_amp = avg_distractor_amp(plot_channel_indices); % [1 x selected_channels]
    plot_no_distractor_amp = avg_no_distractor_amp(plot_channel_indices); % [1 x selected_channels]
    
    plot_chanlocs = cfg.chanlocs(plot_channel_indices);
    
    % Define subplot position
    subplot_idx1 = (i-1)*num_conditions + 1; % Distractor
    subplot_idx2 = (i-1)*num_conditions + 2; % No Distractor
    
    % Plot Distractor Condition
    subplot(num_timepoints, num_conditions, subplot_idx1);
    topoplot(plot_distractor_amp, plot_chanlocs, 'maplimits', fixed_color_limits);
    colorbar;
    title(sprintf('Distractor\nTime = %.2f s', time_points(i)));
    
    % Plot No Distractor Condition
    subplot(num_timepoints, num_conditions, subplot_idx2);
    topoplot(plot_no_distractor_amp, plot_chanlocs, 'maplimits', fixed_color_limits);
    colorbar;
    title(sprintf('No Distractor\nTime = %.2f s', time_points(i)));
end

% Add an overall title to the figure
sgtitle('Topoplots for PO7 and PO8 Across Time Points and Conditions', 'FontSize', 16);

