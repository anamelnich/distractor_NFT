function plotEOG(EOG, labels, params, session, figurepath)
% plotEOGHeatmap plots heat maps of baseline-corrected EOG data.
%
% Inputs:
%   EOG           - 3D matrix [time points x 2 channels x trials]
%   labels        - Vector of trial labels (e.g., 0 for no distractor,
%                   1 and 2 for distractor conditions)

%
% The function performs baseline correction for each trial and channel,
% then creates a heat map (with time on the x-axis and trials on the y-axis)
% for each unique label, plotting horizontal (channel 1) and vertical (channel 2)
% in separate subplots.

    % Baseline correction: for each trial and channel, subtract the mean in baselinePeriod
    baselinePeriod = [-0.1 0];
    baselineIndices = find(params.epochTime >= baselinePeriod(1) & params.epochTime <= baselinePeriod(2));
    EOG_corr = EOG;  % initialize corrected data
    for tr = 1:size(EOG,3)
        for ch = 1:size(EOG,2)
            baseVal = mean(EOG(baselineIndices, ch, tr));
            EOG_corr(:,ch,tr) = EOG(:,ch,tr) - baseVal;
        end
    end
    % Convert amplitudes from microvolts to millivolts
%     EOG_corr = EOG_corr * 1e-3;
    % Get unique labels and loop over them
    plotIndices = find(params.epochTime >= baselinePeriod(1) & params.epochTime <= 0.8);
    plotTime = params.epochTime(plotIndices);
    labels (labels == 2) = 1;
    uniqueLabels = unique(labels);
    for li = 1:length(uniqueLabels)
        currLabel = uniqueLabels(li);
        trialIdx = find(labels == currLabel);

        % Horizontal channel
        horData = squeeze(EOG_corr(:, 2, trialIdx))';  % trials x time
        % Vertical channel
        verData = squeeze(EOG_corr(:, 1, trialIdx))';  % trials x time
        
        % Create a figure for the current label
        figure('Name', sprintf('EOG Heatmap - Label %d', currLabel));
        
        % Plot horizontal channel heat map
        subplot(1,2,1);
        imagesc(plotTime, 1:size(horData,1), horData);
        set(gca, 'YDir', 'normal'); % so trial 1 is at the top
        caxis([-25 25]);
        xline(0);
        xlabel('Time (s)');
        ylabel('Trials');
        title(sprintf('Horizontal EOG (Label %d)', currLabel));
        
        % Plot vertical channel heat map
        subplot(1,2,2);
        imagesc(plotTime, 1:size(verData,1), verData);
        set(gca, 'YDir', 'normal');
        caxis([-25 25]);
        c = colorbar;
        hLabel = ylabel(c, 'Amplitude (\muV)');
        hLabel.Rotation = 270; 
        xline(0);
        xlabel('Time (s)');
        ylabel('Trials');
        title(sprintf('Vertical EOG (Label %d)', currLabel));
        
        safeTitle = regexprep(session, '[^a-zA-Z0-9]', '_');
        label = sprintf('Label %d', currLabel)
        filename = fullfile(figurepath, [safeTitle label ' EOGheat.png']);
        saveas(gcf, filename);
    end
end
