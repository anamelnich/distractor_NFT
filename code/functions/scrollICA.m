function scrollICA(icasig, selectedComponents, datasetName)

    % Basic parameters
    windowSize = 1e4;      % number of samples per window
    startIndexInitial = 2e4;
    
    % Create a figure
    hFig = figure('Name', [datasetName, ' - ICA Components']);

    % Create an axes for plotting
    hAx = axes('Parent', hFig);
    hold(hAx, 'on');
    title(hAx, ['ICA Components for ', datasetName]);
    xlabel(hAx, 'Time (Samples)');
    ylabel(hAx, 'Amplitude (µV) + Offset');

    % Plot first window
    plotWindow(startIndexInitial, windowSize);

    % Create a slider at the bottom of the figure
    hSlider = uicontrol('Parent', hFig, 'Style', 'slider', ...
                        'Units', 'normalized', ...
                        'Position', [0.1 0.01 0.8 0.05], ...
                        'Min', 1, ...
                        'Max', (size(icasig,2) - windowSize), ...
                        'Value', startIndexInitial, ...
                        'Callback', @sliderCallback);

    % Callback for slider movement
    function sliderCallback(hObject, ~)
        newStartIndex = round(hObject.Value);
        cla(hAx); % clear axes
        plotWindow(newStartIndex, windowSize);
    end

    % Helper function to plot a given time window
    function plotWindow(startIdx, winSize)
        timeRange = startIdx : (startIdx + winSize);

        for j = 1:length(selectedComponents)
            compIndex = selectedComponents(j);
            compData = icasig(compIndex, timeRange);

            % Offset each component so it appears on its own row
            % Here, (j-1)*something is the vertical shift between components.
            offsetData = compData + (j - 1) * 10; 
            % ^ '400' is just an example row spacing in microvolts; 
            %   you can adjust for more or less vertical separation.

            plot(hAx, offsetData, 'DisplayName', ['IC' num2str(compIndex)]);
        end

        grid(hAx, 'on');
        
        % For nice labeling along y-axis:
        %
        % If each row is spaced by 400 µV, then:
        % - Row 1 is near 0 µV offset
        % - Row 2 is near 400 µV offset
        % - Row 3 is near 800 µV offset, etc.
        %
        % We can place ticks at these offsets (and label them "IC#"):
        yTickPositions = (0:(length(selectedComponents)-1)) * 10;
        yticks(hAx, yTickPositions);
        yticklabels(hAx, arrayfun(@(x) ['IC' num2str(x)], selectedComponents, 'UniformOutput', false));
        
        % Force x-limits for clarity 
        xlim(hAx, [1, winSize+1]);

        % Force the y-limits so each "row" is limited around ±200 µV 
        % relative to its offset.
        %
        % The bottom row is offset by 0, so we'll set the lower limit to -200 µV
        % The top row is offset by (length(selectedComponents)-1)*400
        % so we add ±50 µV around that.
        topOffset = (length(selectedComponents)-1)*10;
        ylim(hAx, [-0.5, topOffset + 0.5]);
        
    end

end
