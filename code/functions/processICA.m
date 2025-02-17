function dataStruct = processICA(dataStruct)
    % processICA - Perform ICA and visualize components for a single EEG dataset
    %
    % Usage:
    %   processICA(training1)
    %
    %   This function:
    %   1. Runs ICA using 'fastica' on the provided dataset.
    %   2. Plots selected ICA components (1 to 15).
    %   3. Prompts you to specify components to remove.
    %   4. Plots raw signal before and after component removal.
    %   5. Updates the original workspace variable with the cleaned data.
    %
    % Input:
    %   dataStruct - A struct with field 'data' of size [timepoints x channels].
    %                The variable name in the caller workspace is automatically detected.

    %----------------------------------------------------------------------
    % STEP 0: Get the name of the dataset variable from the caller workspace
    %----------------------------------------------------------------------
    datasetName = inputname(1);  % e.g., 'training1'
    if isempty(datasetName)
        error('You must call processICA with a variable, e.g. processICA(training1).');
    end

    disp(['Processing ', datasetName, '...']);

    %----------------------------------------------------------------------
    % STEP 1: Perform ICA
    %----------------------------------------------------------------------
    [icasig, A, W] = fastica(dataStruct.data', ...
        'approach', 'symm', ...
        'g', 'pow3', ...
        'numOfIC', size(dataStruct.data, 2));

    %----------------------------------------------------------------------
    % STEP 2: Plot initial ICA components
    %% ----------------------------------------------------------------------
    selectedComponents = 1:25;
    datasetName = 'EEG';
    scrollICA(icasig, selectedComponents, datasetName);
    % figure('Name', [datasetName, ' - ICA Components']);
    % selectedComponents = 1:15;   % Indices of components to plot
    % windowSize = 1e4;
    % startIndex = 3e4;           % Starting index for the data
    % globalMaxAbsVal = max(max(abs(icasig(selectedComponents, ...
    %                       startIndex : startIndex + windowSize))));
    % 
    % hold on;
    % for j = 1:length(selectedComponents)
    %     compIndex = selectedComponents(j);
    % 
    %     compData = icasig(compIndex, startIndex:startIndex+windowSize);
    % 
    %     % Calculate the max absolute value for THIS component
    %     compMaxVal = max(abs(compData));
    % 
    %     % Normalize each component to +/- 1, then offset by (j-1)
    %     % so each component is plotted in its own "row"
    %     % plot(compData / compMaxVal + (j - 1));
    %     plot(compData / globalMaxAbsVal + (j - 1));
    %     % plot(compData+(j - 1));
    % 
    % end
    % hold off;
    % xlabel('Time (Samples)');
    % ylabel('Amplitude + Offset');
    % title(['ICA Components for ', datasetName]);
    % yticks(0:(length(selectedComponents)-1));
    % yticklabels(arrayfun(@(x) ['IC' num2str(x)], selectedComponents, 'UniformOutput', false));
    % 
    % grid on;

    %% ----------------------------------------------------------------------
    % STEP 3: Manual input for components to remove
    %----------------------------------------------------------------------
    removeIndex = input('Enter the components to remove as a vector (e.g., [8, 10]): ');

    %----------------------------------------------------------------------
    % STEP 4: Plot raw signal before cleaning
    %----------------------------------------------------------------------

    scrollRawEEG(dataStruct, 'preICA');
    % figure('Name', [datasetName, ' - Raw EEG Before ICA']);
    % gap = 50;
    % hold on;
    % plot(bsxfun(@plus, dataStruct.data(3e4:3e4+1e4, :), 0:gap:gap*(size(dataStruct.data, 2)-1)));
    % yLim = get(gca, 'ylim');
    % ylim(yLim);
    % xlabel('Time (Samples)');
    % ylabel('Amplitude (\muV)');
    % title(['Raw EEG Before ICA for ', datasetName]);
    % hold off;

    %----------------------------------------------------------------------
    % STEP 5: Remove specified components and update data
    %----------------------------------------------------------------------
    if ~isempty(removeIndex)
        icasig(removeIndex, :) = 0;              % Zero out the chosen components
        cleanedSignal = A * icasig;              % Reconstruct
        dataStruct.data = cleanedSignal';        % Timepoints x Channels
    end

    %----------------------------------------------------------------------
    % STEP 6: Plot raw signal after cleaning
    %----------------------------------------------------------------------
    scrollRawEEG(dataStruct, 'postICA');
    % figure('Name', [datasetName, ' - Raw EEG After ICA']);
    % hold on;
    % plot(bsxfun(@plus, dataStruct.data(3e4:3e4+1e4, :), 0:gap:gap*(size(dataStruct.data, 2)-1)));
    % yLim = get(gca, 'ylim');
    % ylim(yLim);
    % xlabel('Time (Samples)');
    % ylabel('Amplitude (\muV)');
    % title(['Raw EEG After ICA for ', datasetName]);
    % hold off;

    %----------------------------------------------------------------------
    % STEP 7: Save cleaned data back to the ORIGINAL variable in the base workspace
    %----------------------------------------------------------------------
    % evalin('base', [datasetName ' = dataStruct;']);
    % disp(['Finished processing ', datasetName, '.']);

end
