function combinedEpochs = combineEpochs(epochStructs)
% combineEpochStructs concatenates multiple epoch structures.
%
% Each input epoch structure must have the following fields:
%   - data   : [time x channels x trials]
%   - labels : [trials x 1]
%   - file_id: [trials x 1] (run identifier; e.g., 1 for first run, 2 for second, etc.)
%   - eof    : a vector of end-of-file indices (one per run)
%
% When concatenating, file_id values for each additional epoch are 
% incremented by the maximum file_id from the previous data.
% Similarly, eof values are adjusted by adding the total number of trials
% from the previous epochs.
%
% Input:
%   epochStructs - a cell array containing epoch structures
%
% Output:
%   combinedEpochs - a structure with fields:
%       .data, .labels, .file_id, .eof

combinedEpochs = struct();
combinedEpochs.data = [];       % 3D: time x channels x trials
combinedEpochs.labels = [];     % Column vector, one row per trial
combinedEpochs.file_id = [];    % Column vector of file/run IDs
combinedEpochs.eof = [];        % Vector of end-of-file indices

file_id_offset = 0;


for i = 1:length(epochStructs)
    curEpoch = epochStructs{i};
    
    % Number of trials in current epoch structure:
    nTrials = size(curEpoch.data, 3);
    
    % Concatenate data along the trial (3rd) dimension.
    if isempty(combinedEpochs.data)
        combinedEpochs.data = curEpoch.data;
    else
        combinedEpochs.data = cat(3, combinedEpochs.data, curEpoch.data);
    end
    
    % Concatenate labels (vertically).
    combinedEpochs.labels = [combinedEpochs.labels; curEpoch.labels];
    
    % Adjust file_id values: add the current offset.
    adjusted_file_id = curEpoch.file_id + file_id_offset;
    combinedEpochs.file_id = [combinedEpochs.file_id; adjusted_file_id];
    
    % Update file_id_offset: assume file_id are positive integers; add
    % the maximum file_id from the current structure.
    if ~isempty(curEpoch.file_id)
        file_id_offset = file_id_offset + max(curEpoch.file_id);
    end
    
    %EOF don't need to be adjusted since i only use it for len(eof)
    combinedEpochs.eof = [combinedEpochs.eof; curEpoch.eof];

end

end
