function plotStroopBeh(behSession, plotTitle, figpath)
% plotStroopBeh plots Stroop behavioral data as a bar graph with error bars.
%
% Usage:
%   plotStroopBeh(behSession, plotTitle, figpath)
%
% Inputs:
%   behSession - A structure containing precomputed reaction time statistics:
%                RTn_mean, RTn_std for Neutral,
%                RTc_mean, RTc_std for Congruent,
%                RTi_mean, RTi_std for Incongruent.
%   plotTitle  - A string to use as the figure title and filename.
%   figpath    - The directory where the figure should be saved.
%
% The function plots three bars (Congruent, Neutral, Incongruent) with
% different colors, error bars representing the standard deviations,
% numerical mean values displayed above each bar (offset above the error bars),
% and the Stroop Interference effect (Congruent mean minus Incongruent mean)
% shown at the top. Y-axis limits are set to [400, 1200] ms.

% Reorder data: x = [Congruent, Neutral, Incongruent]
means = [behSession.RTc_mean, behSession.RTn_mean, behSession.RTi_mean];
stds  = [behSession.RTc_std,  behSession.RTn_std,  behSession.RTi_std];

% Define colors: Congruent (blue), Neutral (green), Incongruent (red)
colors = [0, 0.4470, 0.7410;      % Congruent
          0.4660, 0.6740, 0.1880; % Neutral
          0.8500, 0.3250, 0.0980];% Incongruent

x = 1:3;  % x positions for the bars

figure;
hold on;

% Plot the bars with custom colors.
b = bar(x, means, 'FaceColor', 'flat', 'BarWidth', 0.6);
for i = 1:length(x)
    b.CData(i, :) = colors(i, :);
end

% Add error bars (in black)
errorbar(x, means, stds, 'k', 'linestyle', 'none', 'LineWidth', 2);

% Display the numerical mean above each bar.
% Here, the y-position is the mean plus the standard deviation plus an extra offset.
offset = 20;  % extra offset (ms)
for i = 1:length(x)
    yPos = means(i) - stds(i) - offset;
    text(x(i), yPos, sprintf('%.1f', means(i)), ...
         'HorizontalAlignment', 'center', 'FontSize', 12, 'FontWeight', 'bold');
end

% Calculate Stroop Interference effect: Congruent - Incongruent.
interference = means(3) - means(1);
% Display the interference effect at the top center.
text(2, 1150, sprintf('Stroop Interference Effect: %.1f ms', interference), ...
     'HorizontalAlignment', 'center', 'FontSize', 14, 'FontWeight', 'bold', 'Color', [0 0 0]);

% Set axis properties.
set(gca, 'XTick', x, 'XTickLabel', {'Congruent','Neutral','Incongruent'}, 'FontSize', 12);
ylabel('Reaction Time (ms)', 'FontSize', 14);
title(plotTitle, 'Interpreter', 'none', 'FontSize', 16);
ylim([400, 1200]);
grid on;
hold off;

% Save the figure if figpath is provided.
if ~isempty(figpath)
    if ~isfolder(figpath)
        mkdir(figpath);
    end
    % Replace spaces with underscores for the filename.
    fileName = [strrep(plotTitle, ' ', '_') '.png'];
    saveas(gcf, fullfile(figpath, fileName));
end

end


