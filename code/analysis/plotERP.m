function plotERP(data, params, titleName,path)
    % Function to plot ERP data for different conditions and include contra vs none plot
    %
    % Inputs:
    %   data - Struct containing grand average data for different conditions
    %   cfg - Configuration structure with plot parameters
    %   titleName - Title for the figure

    % Plot ERP Distractor Lefts vs Right vs No Distractor
    % figure('Name', [titleName, ' - Left vs Right vs No Distractor']);
    % hold on;
    % plot(params.epochTime, data.grand_average_contra_dright, 'Color', params.plotColor{1}, params.plotOption{:});
    % plot(params.epochTime, data.grand_average_contra_dleft, 'Color', params.plotColor{2}, params.plotOption{:});
    % plot(params.epochTime, data.grand_average_dnone, 'Color', params.plotColor{5}, params.plotOption{:});
    % x_shade = [0.15, 0.3, 0.3, 0.15];
    % y_shade = [-15, -15, 15, 15];
    % %patch(x_shade, y_shade, [0.9 0.9 0.9], 'EdgeColor', 'none', 'FaceAlpha', 0.5);
    % xlabel('Time (s)');
    % ylabel('Amplitude (\muV)');
    % xline(0);
    % yline(0);
    % %xlim([-0.2 0.6]); 
    % ylim([-10 10]);
    % legend('right contra', 'left contra', 'none');
    % title([titleName, ' - Left vs Right vs No Distractor']);

    % Plot ERP Distractor vs No Distractor (Contra vs None)
    figure('Name', [titleName, ' - Distractor vs No Distractor']);
    hold on;
    plot(params.epochTime, data.grand_average_contra, 'Color', params.plotColor{1}, params.plotOption{:});
    plot(params.epochTime, data.grand_average_dnone, 'Color', params.plotColor{5}, params.plotOption{:});
    x_shade = [0.15, 0.5, 0.5, 0.15];
    y_shade = [-15, -15, 15, 15];
    %patch(x_shade, y_shade, [0.9 0.9 0.9], 'EdgeColor', 'none', 'FaceAlpha', 0.5);
    xlabel('Time (s)');
    ylabel('Amplitude (\muV)');
    xline(0);
    yline(0);
    ylim([-10 10]);
    xticks([0 0.3 0.5 1]);
    legend('distractor', 'no distractor');
    title([titleName, ' - Distractor vs No Distractor']);

    safeTitle = regexprep(titleName, '[^a-zA-Z0-9]', '_');  % Remove special characters
    filename = fullfile(path, [safeTitle 'grandavgERP.png']);
    saveas(gcf, filename);
    
    % Plot ERP Distractor Difference wave (Contra - Ipsi)
    figure('Name', [titleName, ' - Distractor vs No Distractor']);
    hold on;
    plot(params.epochTime, data.grand_average_contra, 'Color', params.plotColor{1}, params.plotOption{:});
    plot(params.epochTime, data.grand_average_ipsi, 'Color', params.plotColor{2}, params.plotOption{:});
    plot(params.epochTime, (data.grand_average_contra - data.grand_average_ipsi), 'Color', params.plotColor{6}, params.plotOption{:});
    x_shade = [0.15, 0.5, 0.5, 0.15];
    y_shade = [-15, -15, 15, 15];
    %patch(x_shade, y_shade, [0.9 0.9 0.9], 'EdgeColor', 'none', 'FaceAlpha', 0.5);
    xlabel('Time (s)');
    ylabel('Amplitude (\muV)');
    xline(0);
    yline(0);
    ylim([-10 10]);
    xticks([0 0.3 0.5 1]);
    legend('contra', 'ipsi', 'difference');
    title([titleName, ' - Contra - Ipsi Difference Wave']);

    safeTitle = regexprep(titleName, '[^a-zA-Z0-9]', '_');  % Remove special characters
    filename = fullfile(path, [safeTitle 'DifferenceWaveD.png']);
    saveas(gcf, filename);

    
    % Plot ERP No Distractor Difference wave
    figure('Name', [titleName, ' - Distractor vs No Distractor']);
    hold on;
    plot(params.epochTime, data.grand_average_dnone_right, 'Color', params.plotColor{1}, params.plotOption{:});
    plot(params.epochTime, data.grand_average_dnone_left, 'Color', params.plotColor{2}, params.plotOption{:});
    plot(params.epochTime, (data.grand_average_dnone_right - data.grand_average_dnone_left), 'Color', params.plotColor{6}, params.plotOption{:});
    plot(params.epochTime, (data.grand_average_dnone_diff), 'Color', params.plotColor{5}, params.plotOption{:});
    x_shade = [0.15, 0.5, 0.5, 0.15];
    y_shade = [-15, -15, 15, 15];
    %patch(x_shade, y_shade, [0.9 0.9 0.9], 'EdgeColor', 'none', 'FaceAlpha', 0.5);
    xlabel('Time (s)');
    ylabel('Amplitude (\muV)');
    xline(0);
    yline(0);
    ylim([-10 10]);
    xticks([0 0.3 0.5 1]);
    legend('right electrodes', 'left electrodes', 'right - left', 'random side diff');
    title([titleName, ' - No Distractor Difference Wave']);
    safeTitle = regexprep(titleName, '[^a-zA-Z0-9]', '_');  % Remove special characters
    filename = fullfile(path, [safeTitle 'DifferenceWaveND.png']);
    saveas(gcf, filename);
    
    % Plot ERP Distractor vs No Distractor Difference wave
    figure('Name', [titleName, ' - Distractor vs No Distractor Difference Waves']);
    hold on;
    plot(params.epochTime, data.grand_average_d_diff, 'Color', params.plotColor{1}, params.plotOption{:});
    plot(params.epochTime, data.grand_average_dnone_diff, 'Color', params.plotColor{5}, params.plotOption{:});
    x_shade = [0.15, 0.5, 0.5, 0.15];
    y_shade = [-15, -15, 15, 15];
    %patch(x_shade, y_shade, [0.9 0.9 0.9], 'EdgeColor', 'none', 'FaceAlpha', 0.5);
    xlabel('Time (s)');
    ylabel('Amplitude (\muV)');
    xline(0);
    yline(0);
    ylim([-6 6]);
    xTickInterval = 0.1;  % Tick marks every 1 second
    xticks(0:xTickInterval:max(params.epochTime));
    legend('Distractor contra-ipsi', 'No Distractor random side subtraction');
    title([titleName, ' - Distractor vs No distratcor Difference Wave']);
    safeTitle = regexprep(titleName, '[^a-zA-Z0-9]', '_');  % Remove special characters
    filename = fullfile(path, [safeTitle 'DifferenceWaveDvsND.png']);
    saveas(gcf, filename);

    
end