function plotBehavior(beh, cfg, titleName, path)
    % Plot RT distributions
    figure;
    hold on;

    hRTnd = histogram(beh.RTnd, 'FaceColor', cfg.plotColor{5}, 'FaceAlpha', 0.5,'BinWidth', 40);
    hRTd = histogram(beh.RTd, 'FaceColor', cfg.plotColor{1}, 'FaceAlpha', 0.5,'BinWidth', 40);
        
    xline(beh.RTnd_mean, 'Color', cfg.plotColor{5}, 'LineWidth', 2);
    xline(beh.RTd_mean, 'Color', cfg.plotColor{1}, 'LineWidth', 2);

    xlabel('Reaction Time (ms)');
    ylabel('Counts');
    xlim([200 1400]);
    ylim([0 40]);
    
    leg1 = sprintf('No Distractor: %.0f ± %.0f ms', beh.RTnd_mean, beh.RTnd_std);
    leg2 = sprintf('Distractor: %.0f ± %.0f ms', beh.RTd_mean, beh.RTd_std);
    legend([hRTnd, hRTd], {leg1, leg2});
    
    title([titleName, ' - Reaction Times']);
    safeTitle = regexprep(titleName, '[^a-zA-Z0-9]', '_');  % Remove special characters
    filename = fullfile(path, [safeTitle 'RT.png']);
    saveas(gcf, filename);
    hold off;
    
     % Plot RT distributions for correct responses only 
    figure;
    hold on;

    hRTnd = histogram(beh.RTnd_correct, 'FaceColor', cfg.plotColor{5}, 'FaceAlpha', 0.5,'BinWidth', 40);
    hRTd = histogram(beh.RTd_correct, 'FaceColor', cfg.plotColor{1}, 'FaceAlpha', 0.5,'BinWidth', 40);
        
    xline(beh.RTnd_mean_corr, 'Color', cfg.plotColor{5}, 'LineWidth', 2);
    xline(beh.RTd_mean_corr, 'Color', cfg.plotColor{1}, 'LineWidth', 2);
    
    xlabel('Reaction Time (ms)');
    ylabel('Counts');
    xlim([200 1400]);
    ylim([0 40]);
    
    leg1 = sprintf('No Distractor: %.0f ± %.0f ms', beh.RTnd_mean_corr, beh.RTnd_std_corr);
    leg2 = sprintf('Distractor: %.0f ± %.0f ms', beh.RTd_mean_corr, beh.RTd_std_corr);
    legend([hRTnd, hRTd], {leg1, leg2});
 
    title([titleName, ' - Reaction Times (correct only)']);
    safeTitle = regexprep(titleName, '[^a-zA-Z0-9]', '_');  % Remove special characters
    filename = fullfile(path, [safeTitle 'RTcorr.png']);
    saveas(gcf, filename);
    hold off;
end