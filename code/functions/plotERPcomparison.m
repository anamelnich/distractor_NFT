function plotERPcomparison(processed_data, params, titleName, path, sessions)
    % Number of sessions in the list
    numSessions = length(sessions);
    % Compute the number of complete pairs (if odd number, the last one is ignored)
    nPairs = floor(numSessions / 2);
    
    %% Loop over pairs for ERP Grand Averages
    for p = 1:nPairs
        % Determine the indices for the current pair
        idx1 = 2*p - 1;
        idx2 = 2*p;
        
        % Create a new figure for this pair
        figure('Name', [titleName, ' - Pair ', num2str(p), ' Distractor vs No Distractor']);
        hold on;
        
        % Add shading for a specific time window (optional)
        x_shade = [0.15, 0.5, 0.5, 0.15];
        y_shade = [-15, -15, 15, 15];
        patch(x_shade, y_shade, [0.9 0.9 0.9], 'EdgeColor', 'none', ...
            'FaceAlpha', 0.5, 'HandleVisibility', 'off');
        
        legendEntries = {};
        
        % Loop through the two sessions in the current pair
        for j = [idx1, idx2]
            sessName = sessions{j};
            if isfield(processed_data, sessName)
                sessData = processed_data.(sessName).ROIerp;
                
                % For the first session in the pair use special colors,
                % for the second session use default colors.
                if j == idx1
                    color_contra = params.plotColor{3};
                    color_dnone  = params.plotColor{7};
                else
                    color_contra = params.plotColor{1};
                    color_dnone  = params.plotColor{5};
                end
                
                lineWidth = 2;  % Adjust line width as needed
                
                % Plot grand average for distractor (contra)
                plot(params.epochTime, sessData.grand_average_contra, ...
                    'Color', color_contra, 'LineWidth', lineWidth);
                legendEntries{end+1} = [sessName, ' distractor'];
                
                % Plot grand average for no distractor
                plot(params.epochTime, sessData.grand_average_dnone, ...
                    'Color', color_dnone, 'LineWidth', lineWidth);
                legendEntries{end+1} = [sessName, ' no distractor'];
            end
        end
        
        xlabel('Time (s)');
        ylabel('Amplitude (\muV)');
        xline(0, '--k');
        yline(0, '--k');
        ylim([-6 10]);
        xlim([-0.1 0.6]);
        xticks(0:0.1:max(params.epochTime));
        legend(legendEntries, 'Location', 'best');
        title([titleName, ' - Pair ', num2str(p), ' Distractor vs No Distractor']);
        
        % Save the figure with pair indicated in the filename
        safeTitle = regexprep([titleName, '_Pair_', num2str(p)], '[^a-zA-Z0-9]', '_');
        filename = fullfile(path, [safeTitle, 'grandavgERP.png']);
        saveas(gcf, filename);
    end

    %% Loop over pairs for ERP Difference Waves
    for p = 1:nPairs
        idx1 = 2*p - 1;
        idx2 = 2*p;
        
        figure('Name', [titleName, ' - Pair ', num2str(p), ' Distractor vs No Distractor Difference Waves']);
        hold on;
        
        x_shade = [0.15, 0.5, 0.5, 0.15];
        y_shade = [-15, -15, 15, 15];
        patch(x_shade, y_shade, [0.9 0.9 0.9], 'EdgeColor', 'none', ...
            'FaceAlpha', 0.5, 'HandleVisibility', 'off');
        
        legendEntries = {};
        for j = [idx1, idx2]
            sessName = sessions{j};
            if isfield(processed_data, sessName)
                sessData = processed_data.(sessName).ROIerp;
                
                if j == idx1
                    color_contra = params.plotColor{3};
                    color_dnone  = params.plotColor{7};
                else
                    color_contra = params.plotColor{1};
                    color_dnone  = params.plotColor{5};
                end
                
                lineWidth = 2;
                
                % Plot difference wave for distractor
                plot(params.epochTime, sessData.grand_average_d_diff, ...
                    'Color', color_contra, 'LineWidth', lineWidth);
                legendEntries{end+1} = [sessName, ' distractor'];
                
                % Plot difference wave for no distractor
                plot(params.epochTime, sessData.grand_average_dnone_diff, ...
                    'Color', color_dnone, 'LineWidth', lineWidth);
                legendEntries{end+1} = [sessName, ' no distractor'];
            end
        end
        
        xlabel('Time (s)');
        ylabel('Amplitude (\muV)');
        xline(0, '--k');
        yline(0, '--k');
        ylim([-4 4]);
        xlim([-0.1 0.6]);
        xticks(0:0.1:max(params.epochTime));
        legend(legendEntries, 'Location', 'best');
        title([titleName, ' - Pair ', num2str(p), ' Distractor vs No Distractor Difference Waves']);
        
        safeTitle = regexprep([titleName, '_Pair_', num2str(p)], '[^a-zA-Z0-9]', '_');
        filename = fullfile(path, [safeTitle, 'DifferenceWaveDvsND.png']);
        saveas(gcf, filename);
    end
end

