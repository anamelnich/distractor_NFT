function conversionFactors = convertEOGtoAngles(EOG, labels, cfg, session, figpath)
% convertEOGtoAngles: Convert bipolar EOG signals (in microvolts) to eye movement
% angles (in visual degrees) and create several diagnostic plots.
%
% This function:
%   1. Performs baseline correction (using -0.1 to 0 s) on the EOG data.
%   2. Computes conversion factors (microvolts/degree) separately for vertical and
%      horizontal channels based on trials where the target was at 5° visual angle.
%      For vertical: label 1 (top, +5°) and label 3 (bottom, -5°).
%      For horizontal: label 2 (right, +5°) and label 4 (left, -5°).
%      A linear regression forced through the origin is used.
%   3. Creates a set of figures:
%       a. For each unique label, a figure is generated with two subplots
%          (vertical and horizontal EOG) showing individual trial waveforms and
%          an overlaid average waveform.
%       b. A 2D trajectory plot is generated in which at each time point the
%          horizontal and vertical positions (converted to degrees) are plotted.
%
% INPUTS:
%   EOG      - 3D matrix [time points x 2 channels x trials] (in microvolts)
%   labels   - Vector of trial labels (1: top, 2: right, 3: bottom, 4: left)
%   cfg      - Configuration structure containing a field 'epochTime' (time vector)
%   session  - String identifier for the session (e.g., [subjectID ' ' fname])
%   figpath  - Folder path where figures will be saved
%
% OUTPUT:
%   conversionFactors - Structure with fields:
%       .vertical   - microvolts per visual degree (vertical channel)
%       .horizontal - microvolts per visual degree (horizontal channel)
%


%% Baseline correction
baselinePeriod = [-0.1 0];
baselineIdx = find(cfg.epochTime >= baselinePeriod(1) & cfg.epochTime <= baselinePeriod(2));
EOG_corr = EOG;  % initialize corrected data
nTrials = size(EOG, 3);
for tr = 1:nTrials
    for ch = 1:size(EOG,2)
        baseVal = mean(EOG(baselineIdx, ch, tr));
        EOG_corr(:, ch, tr) = EOG(:, ch, tr) - baseVal;
    end
end

%% Define response period for analysis and get plotting time vector
respPeriod = [0 0.4];  % adjust as needed
respIdx = find(cfg.epochTime >= respPeriod(1) & cfg.epochTime <= respPeriod(2));
plotTime = cfg.epochTime(respIdx);

%% Regression: Compute conversion factors (microvolts per degree)
% For vertical channel (channel 1): Labels 1 (top, +5°) and 3 (bottom, -5°)
vertical_trials = find(labels == 1 | labels == 3);
vertical_amp = zeros(length(vertical_trials),1);
vertical_target = zeros(length(vertical_trials),1);
for i = 1:length(vertical_trials)
    tr = vertical_trials(i);
    if labels(tr) == 1  % top: expect positive deflection
        vertical_amp(i) = max(EOG_corr(respIdx, 1, tr));
        vertical_target(i) = 5;  % degrees
    elseif labels(tr) == 3  % bottom: expect negative deflection
        vertical_amp(i) = min(EOG_corr(respIdx, 1, tr));
        vertical_target(i) = -5;
    end
end
% Linear regression forced through the origin: slope = sum(amp*target)/sum(target^2)
slope_vertical = (vertical_amp' * vertical_target) / (vertical_target' * vertical_target);

% For horizontal channel (channel 2): Labels 2 (right, +5°) and 4 (left, -5°)
horizontal_trials = find(labels == 2 | labels == 4);
horizontal_amp = zeros(length(horizontal_trials),1);
horizontal_target = zeros(length(horizontal_trials),1);
for i = 1:length(horizontal_trials)
    tr = horizontal_trials(i);
    if labels(tr) == 2  % right: expect positive deflection
        horizontal_amp(i) = max(EOG_corr(respIdx, 2, tr));
        horizontal_target(i) = 5;
    elseif labels(tr) == 4  % left: expect negative deflection
        horizontal_amp(i) = min(EOG_corr(respIdx, 2, tr));
        horizontal_target(i) = -5;
    end
end
slope_horizontal = (horizontal_amp' * horizontal_target) / (horizontal_target' * horizontal_target);

% Take absolute values (conversion factor: microvolts per degree)
convFactor.vertical = abs(slope_vertical);
convFactor.horizontal = abs(slope_horizontal);
conversionFactors = convFactor;  % output structure

fprintf('Vertical conversion factor: %.2f microvolts/degree\n', convFactor.vertical);
fprintf('Horizontal conversion factor: %.2f microvolts/degree\n', convFactor.horizontal);

%% Plot individual trial waveforms for each label
uniqueLabels = unique(labels);
for li = 1:length(uniqueLabels)
    currLabel = uniqueLabels(li);
    trialIdx = find(labels == currLabel);
    
    % Extract data for the response period (for each channel)
    vertical_data = squeeze(EOG_corr(respIdx, 1, trialIdx));    % [time x trials]
    horizontal_data = squeeze(EOG_corr(respIdx, 2, trialIdx));  % [time x trials]
    
    % Create figure for current label
    figure('Name', sprintf('EOG Waveforms - Label %d', currLabel));
    
    % Vertical channel subplot
    subplot(2,1,1);
    plot(plotTime, vertical_data, 'Color', [0.7 0.7 0.7]); hold on;
    avgVertical = mean(vertical_data, 2);
    plot(plotTime, avgVertical, 'k', 'LineWidth', 2);
    xlabel('Time (s)'); ylabel('Vertical EOG (\muV)');
    title(sprintf('Label %d - Vertical Channel', currLabel));
    hold off;
    
    % Horizontal channel subplot
    subplot(2,1,2);
    plot(plotTime, horizontal_data, 'Color', [0.7 0.7 0.7]); hold on;
    avgHorizontal = mean(horizontal_data, 2);
    plot(plotTime, avgHorizontal, 'k', 'LineWidth', 2);
    xlabel('Time (s)'); ylabel('Horizontal EOG (\muV)');
    title(sprintf('Label %d - Horizontal Channel', currLabel));
    hold off;
    
    % Save the figure using a safe file name
    safeSession = regexprep(session, '[^a-zA-Z0-9]', '_');
    filename = fullfile(figpath, sprintf('%s_Label%d_EOG_waveforms.png', safeSession, currLabel));
    saveas(gcf, filename);
end

%% Plot 2D trajectory of eye movements (converted to visual degrees)
% % Compute unique labels and prepare colors
% uniqueLabels = unique(labels);
% colors = lines(length(uniqueLabels));
% 
% % Precompute average trajectories for each label (in degrees)
% avgHor_all = cell(length(uniqueLabels),1);
% avgVert_all = cell(length(uniqueLabels),1);
% for li = 1:length(uniqueLabels)
%     currLabel = uniqueLabels(li);
%     trialIdx = (labels == currLabel);
%     % Average over the third dimension (trials) and convert to degrees
%     avgHor_all{li} = mean(EOG_corr(respIdx, 2, trialIdx), 3) / convFactor.horizontal;
%     avgVert_all{li} = mean(EOG_corr(respIdx, 1, trialIdx), 3) / convFactor.vertical;
% end
% 
% % Create figure and subplots arranged horizontally (1 row, multiple columns)
% fig = figure('Name','Animated Average Eye Movement Trajectories');
% animLines = gobjects(length(uniqueLabels),1);
% for li = 1:length(uniqueLabels)
%     subplot(1, length(uniqueLabels), li);  % Arrange subplots horizontally
%     hold on;
%     % Set fixed axis limits to ±8 degrees and force square axes
%     axis([-8 8 -8 8]);
%     axis square;
%     xlabel('Horizontal Angle (°)');
%     ylabel('Vertical Angle (°)');
%     % title(sprintf('Label %d', uniqueLabels(li)));
% 
%     % Create an animated line for the average trajectory
%     animLines(li) = animatedline('Color', colors(li,:), 'LineWidth', 2);
% 
%     % Add fixation cross at (0,0)
%     fix_size = 0.3; % half-length of cross arms
%     plot([-fix_size, fix_size], [0, 0], 'k', 'LineWidth', 1);
%     plot([0, 0], [-fix_size, fix_size], 'k', 'LineWidth', 1);
% 
%     % Add target dot based on label
%     switch uniqueLabels(li)
%         case 1
%             targetX = 0; targetY = 5;
%         case 2
%             targetX = 5; targetY = 0;
%         case 3
%             targetX = 0; targetY = -5;
%         case 4
%             targetX = -5; targetY = 0;
%         otherwise
%             targetX = 0; targetY = 0;
%     end
%     plot(targetX, targetY, 'ko', 'MarkerFaceColor', 'k', 'MarkerSize', 5);
% 
%     hold off;
% end
% 
% % Determine number of time points in the response period
% numPoints = length(respIdx);
% nLoops = 50;  % number of times to loop the animation
% 
% % Animate the average trajectories in each subplot concurrently
% for loopIter = 1:nLoops
%     % Clear previous animated points for each subplot
%     for li = 1:length(uniqueLabels)
%         subplot(1, length(uniqueLabels), li);
%         clearpoints(animLines(li));
%     end
% 
%     % Animate frame-by-frame
%     for j = 1:numPoints
%         for li = 1:length(uniqueLabels)
%             % Add a point from the precomputed average trajectory
%             addpoints(animLines(li), avgHor_all{li}(j), avgVert_all{li}(j));
%         end
%         drawnow;
%     end
%     pause(1);  % pause between animation loops
% end
% 

% Assume uniqueLabels, colors, respIdx, EOG_corr, and convFactor are defined
% Create figure and one set of axes
fig = figure('Name','Animated Trial-by-Trial Eye Movement Trajectories');
hold on;
axis([-8 8 -8 8]);
axis square;
xlabel('Horizontal Angle (°)');
ylabel('Vertical Angle (°)');
title('Animated Trial-by-Trial Eye Movement Trajectories');

% Draw fixation cross at (0,0)
fix_size = 0.3; % half-length of cross arms
plot([-fix_size, fix_size], [0, 0], 'k', 'LineWidth', 1);
plot([0, 0], [-fix_size, fix_size], 'k', 'LineWidth', 1);

% Define unique labels and assign colors based on label value (assumes labels are 1 to 4)
uniqueLabels = unique(labels);
colors = lines(length(uniqueLabels));  % colors(1,:) for label 1, etc.

% Draw target dots for each label
for li = 1:length(uniqueLabels)
    switch uniqueLabels(li)
        case 1, targetX = 0; targetY = 5;
        case 2, targetX = 5; targetY = 0;
        case 3, targetX = 0; targetY = -5;
        case 4, targetX = -5; targetY = 0;
        otherwise, targetX = 0; targetY = 0;
    end
    plot(targetX, targetY, 'ko', 'MarkerFaceColor', 'k', 'MarkerSize', 5);
end

% Number of animation loops (repeat entire set of trials if desired)
nLoops = 3;
numTrials = size(EOG_corr, 3);

for loopIter = 1:nLoops
    % Loop through trials in their original order
    for t = 1:numTrials
        % Get current trial's label (assumed to be 1,2,3,4)
        currLabel = labels(t);
        % Use the label as an index to assign color (if labels are 1-4)
        trialColor = colors(currLabel, :);
        
        % Create an animated line for the current trial
        animLine = animatedline('Color', trialColor, 'LineWidth', 2);
        
        % Animate point-by-point over the response period
        for j = 1:length(respIdx)
            x = EOG_corr(respIdx(j), 2, t) / convFactor.horizontal;
            y = EOG_corr(respIdx(j), 1, t) / convFactor.vertical;
            addpoints(animLine, x, y);
            drawnow;
            pause(0.01);  % adjust pause to control the speed of the animation
        end
        
        pause(0.5);   % pause to view the completed trial trajectory
        delete(animLine);  % erase the line before showing the next trial
    end
end

% Save the 2D trajectory figure
safeSession = regexprep(session, '[^a-zA-Z0-9]', '_');
filename = fullfile(figpath, sprintf('%s_EOG_trajectories.png', safeSession));
saveas(gcf, filename);

end
