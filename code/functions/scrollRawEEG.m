function scrollRawEEG(dataStruct, datasetName)
% SCROLLRAweEG Scroll through raw EEG data in an interactive figure.
%
%   SCROLLRAweEG(dataStruct, datasetName) displays the EEG data in
%   dataStruct.data (assumed to be samples x channels) in an interactive
%   figure with a slider. 
%
%   - dataStruct.data: a [samples x channels] matrix of raw EEG
%   - datasetName:     a string for labeling the figure

    % -- Parameters you can adjust --
    windowSize       = 10e3;   % how many samples to display at a time
    startIndexInitial= 30e3;   % starting sample index for the first window
    gap              = 50;     % vertical offset between channels

    % -- Basic checks --
    [nSamples, nChans] = size(dataStruct.data);
    if startIndexInitial + windowSize > nSamples
        % Clip start index so we don't exceed total data length
        startIndexInitial = max(1, nSamples - windowSize);
    end

    % -- Create a figure and axes --
    hFig = figure('Name',[datasetName,' - Raw EEG Scroller'], ...
                  'Color','white');
    hAx  = axes('Parent',hFig);
    hold(hAx, 'on');
    title(hAx, ['Raw EEG for ', datasetName]);
    xlabel(hAx, 'Time (Samples)');
    ylabel(hAx, 'Amplitude (\muV) + offset');
    set(hAx, 'YDir','reverse'); 
    % Some people prefer 'reverse' so channel 1 is at the top, but it's optional

    % -- Initial plot of the first window --
    plotWindow(startIndexInitial, windowSize);

    % -- Create a slider for scrolling --
    hSlider = uicontrol('Parent',hFig, 'Style','slider', ...
        'Units','normalized', 'Position',[0.1 0.01 0.8 0.05], ...
        'Min', 1, 'Max', nSamples - windowSize, ...
        'Value', startIndexInitial, ...
        'Callback', @sliderCallback);

    % Callback function for the slider
    function sliderCallback(hObject, ~)
        newStartIndex = round(hObject.Value);
        cla(hAx); % clear axes
        plotWindow(newStartIndex, windowSize);
    end

    % Nested function to plot a given window
    function plotWindow(startIdx, winSize)
        endIdx = startIdx + winSize - 1;
        if endIdx > nSamples
            endIdx = nSamples;
        end
        timeRange = startIdx:endIdx;

        % Extract the portion of data you want to plot
        subData = dataStruct.data(timeRange, :); % [winSize x nChans]

        % Offset each channel so they don't overlap
        %   channel #1 has offset 0
        %   channel #2 has offset gap
        %   channel #3 has offset 2*gap
        offsets = 0:gap:(gap*(nChans-1));
        % Add the offsets column-wise
        offsetData = bsxfun(@plus, subData, offsets);

        % Plot each column as a separate line (MATLAB auto interprets row=sample)
        plot(hAx, offsetData, 'LineWidth', 1);
        grid(hAx, 'on');

        % Set x-limits to show exactly [1..windowSize] samples
        xlim(hAx, [1 size(offsetData,1)]);

        % A simple Y-lim that fits all channels; you can also let MATLAB autoscale
        ylims = [min(offsetData(:)) max(offsetData(:))];
        ylim(hAx, ylims);
    end

end
