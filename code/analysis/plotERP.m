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
    x_shade = [0.15, 0.3, 0.3, 0.15];
    y_shade = [-15, -15, 15, 15];
    %patch(x_shade, y_shade, [0.9 0.9 0.9], 'EdgeColor', 'none', 'FaceAlpha', 0.5);
    xlabel('Time (s)');
    ylabel('Amplitude (\muV)');
    xline(0);
    yline(0);
    ylim([-10 10]);
    legend('distractor', 'no distractor');
    title([titleName, ' - Distractor vs No Distractor']);

    safeTitle = regexprep(titleName, '[^a-zA-Z0-9]', '_');  % Remove special characters
    filename = fullfile(path, [safeTitle 'grandavgERP.png']);
    saveas(gcf, filename);
end