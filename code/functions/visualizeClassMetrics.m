function visualizeClassMetrics(behData, params, path)
    % Define decoding sessions
    decodingSessions = {'decoding1', 'decoding2', 'decoding3'};
    nSessions = numel(decodingSessions);

    plotColor = {
    [0.294, 0, 0.510],      % Blue1 (Indigo)
    [0, 0.502, 0],          % Green1 (Dark Green)
    [0.698, 0.133, 0.133],  % Red1 (Dark Red)
    [0, 0, 0.502],          % Blue2 (Navy blue)
    [0.133, 0.545, 0.133],  % Green2 (Forest green)
    [0.698, 0.133, 0.133],  % Red2 (Firebrick)
    [0.255, 0.412, 0.882],  % Blue3 (Royal blue)
    [0.196, 0.804, 0.196],  % Green3 (Lime green)
    [0.804, 0.361, 0.361]   % Red3 (Indian red)
                        };
    
    %% 2. Run-by-Run Metrics Plot (accuracy, TPR, TNR)
    % We'll concatenate runs from each session along the x-axis.
    x_all = [];
    acc_all = [];
    TPR_all = [];
    TNR_all = [];
    sessionBoundaries = [];  % to mark separation between sessions
    offset = 0;
    figure('Name','Run-by-Run Classification Metrics');
    hold on;
    for s = 1:nSessions
        sessionName = decodingSessions{s};
        if isfield(behData, sessionName)
            runMetrics = behData.(sessionName).classMetrics.by_run;
            nRuns = length(runMetrics.runNumbers);
            % x positions for these runs
            x = (1:nRuns) + offset;
            x_all = [x_all; x(:)];
            acc_all = [acc_all; runMetrics.accuracy(:)];
            TPR_all = [TPR_all; runMetrics.TPR(:)];
            TNR_all = [TNR_all; runMetrics.TNR(:)];
            
            acc_color = plotColor{1};
            tpr_color = plotColor{2};
            tnr_color = plotColor{3};

            plot(x, runMetrics.accuracy, '-o', 'LineWidth', 2, 'Color', acc_color);
            plot(x, runMetrics.TPR, '--s', 'LineWidth', 1, 'Color', tpr_color);
            plot(x, runMetrics.TNR, '--d', 'LineWidth', 1, 'Color', tnr_color);
            
            % Add overall accuracy annotation for this session.
            overallAcc = behData.(sessionName).classMetrics.overall.accuracy;
            x_center = offset + nRuns/2;
            % Set y position slightly above the maximum run accuracy for this session.
            allMetrics = [runMetrics.accuracy(:); runMetrics.TPR(:); runMetrics.TNR(:)];
            overallMin = min(allMetrics);
            y_pos = max(overallMin) - 0.02;
            text(x_center, y_pos, sprintf('Session Acc: %.2f', overallAcc), ...
                'Color', acc_color, 'FontSize', 10, 'FontWeight', 'bold', ...
                'HorizontalAlignment', 'center');
            
            offset = offset + nRuns;
            sessionBoundaries(end+1) = offset;
        end
    end
    hold off;
    xlabel('Run');
    ylabel('Percent(%)');
    xTickInterval = 1;  % Tick marks every 1 second
    xticks(0:xTickInterval:23);
    title('Run-by-Run Accuracy, TPR, and TNR');


% Now create the legend using only these dummy handles.
legend({'Accuracy', 'TPR', 'TNR'}, 'Location', 'best');
    % Optionally, mark session boundaries with vertical lines.
    for b = 1:length(sessionBoundaries)-1
        xline(sessionBoundaries(b)+0.5, '--k');
    end
    saveas(gcf, fullfile(path, 'RunByRunMetrics.png'));

    
    %% 3. Confusion Matrices for Each Decoding Session (subplots)
    confusionMats = cell(nSessions, 1);
    for s = 1:nSessions
        sessionName = decodingSessions{s};
        if isfield(behData, sessionName)
            trueLabels = behData.(sessionName).trial_type;
            predLabels = behData.(sessionName).class;
            CM = confusionmat(trueLabels, predLabels);
            confusionMats{s} = CM;
        else
            confusionMats{s} = [];
        end
    end
    
    figure('Name','Confusion Matrices for Decoding Sessions');
    for s = 1:nSessions
        subplot(1, nSessions, s);
        if ~isempty(confusionMats{s})
            CM = confusionMats{s};
            % Compute percentages per row (for each true class)
            CM_percent = zeros(size(CM));
            for i = 1:size(CM,1)
                
                rowSum = sum(CM(i,:));
                if rowSum > 0
                    CM_percent(i,:) = 100 * CM(i,:) / rowSum;
                else
                    CM_percent(i,:) = 0;
                end
            end
    
            imagesc(CM_percent);
            colormap('summer');  
            caxis([30 75]);  % percentages range 0 to 100
            colorbar;
            title(decodingSessions{s});
            xlabel('Predicted Class');
            ylabel('True Class');
            xticks([1 2]);
            xticklabels({'0','1'});
            yticks([1 2]);
            yticklabels({'0','1'});
            axis square;
    
            % Annotate each cell with its percentage value
            [nRows, nCols] = size(CM_percent);
            for i = 1:nRows
                for j = 1:nCols
                    if CM_percent(i,j) > 50
                        txtColor = 'k';
                    else
                        txtColor = 'w';
                    end
                    text(j, i, sprintf('%.1f%%', CM_percent(i,j)), ...
                        'HorizontalAlignment', 'center', 'Color', txtColor, 'FontWeight', 'bold');
                end
            end
        else
            title([decodingSessions{s}, ' - No Data']);
        end
    end

    saveas(gcf, fullfile(path, 'ConfusionMatrices.png'));
end
