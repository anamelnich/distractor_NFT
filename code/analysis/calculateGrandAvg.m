function averages = calculateGrandAvg(data, labels, leftElectrodeIndices, rightElectrodeIndices)
    % Calculate means and grand averages for given epochs and electrode indices.
    %
    % Inputs:
    %   epochs - Struct containing EEG data and labels
    %   leftElectrodeIndices - Array of indices for left electrodes
    %   rightElectrodeIndices - Array of indices for right electrodes
    %
    % Outputs:
    %   results - Struct containing grand averages and intermediate calculations
    
    % Compute data averages for each condition
    leftElect_dright = mean(data(:, leftElectrodeIndices, labels == 1), 3);
    rightElect_dleft = mean(data(:, rightElectrodeIndices, labels == 2), 3);
    leftElect_dnone = mean(data(:, leftElectrodeIndices, labels == 0), 3);
    rightElect_dnone = mean(data(:, rightElectrodeIndices, labels == 0), 3);


    % Grand averages for different conditions
    grand_average_contra_dleft = mean(rightElect_dleft, 2);
    grand_average_contra_dright = mean(leftElect_dright, 2);
    grand_average_contra = mean([rightElect_dleft, leftElect_dright], 2);
    grand_average_dnone_left = mean(leftElect_dnone, 2);
    grand_average_dnone_right = mean(rightElect_dnone, 2);
    grand_average_dnone = mean([rightElect_dnone, leftElect_dnone], 2);

    % Store results in a struct
    averages.leftElect_dright = leftElect_dright;
    averages.rightElect_dleft = rightElect_dleft;
    averages.leftElect_dnone = leftElect_dnone;
    averages.rightElect_dnone = rightElect_dnone;
    averages.grand_average_contra_dleft = grand_average_contra_dleft;
    averages.grand_average_contra_dright = grand_average_contra_dright;
    averages.grand_average_contra = grand_average_contra;
    averages.grand_average_dnone_left = grand_average_dnone_left;
    averages.grand_average_dnone_right = grand_average_dnone_right;
    averages.grand_average_dnone = grand_average_dnone;
end
