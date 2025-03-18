function data = loadData(dataPath, subjectID)
% loadDataScalable loads EEG recordings from multiple day folders.
%
% Assumptions:
%  - Day folders are named as: subjectID_YYYYMMDD
%    e.g., 'e5_20240316' for recordings on March 16, 2024.
%  - Within each day folder, subfolders are named as:
%    subjectID_YYYYMMDDHHMMSS_task
%    e.g., 'e5_20240316114601_training', 'e5_20240316120000_decoding'
%
% Recordings for the same task on the same day are concatenated. For
% different days, the recordings are stored in separate fields.
%
% Example output: data.decoding1, data.decoding2, etc.

% List day folders for the subject (e.g. 'e5_20*')
dayFolders = dir(fullfile(dataPath, [subjectID '_20*']));
if isempty(dayFolders)
    error('No day folders found for subject %s in %s', subjectID, dataPath);
end

% We'll group recordings by task and by day.
% Use a nested structure: tasks.(taskType).(dayField) = cell array of sessions.
tasks = struct();

for d = 1:length(dayFolders)
    dayFolderName = dayFolders(d).name;  % e.g., 'e5_20240316'
    dayFolderPath = fullfile(dayFolders(d).folder, dayFolderName);
    
    % Extract the date from the day folder name.
    % Expected format: subjectID_YYYYMMDD, e.g. 'e5_20240316'
    dayTokens = regexp(dayFolderName, ['^' subjectID '_(\d{8})$'], 'tokens');
    if isempty(dayTokens)
        warning('Day folder name "%s" does not match expected format.', dayFolderName);
        continue;
    end
    % Make a valid field name by prepending a letter (e.g., 'd')
    dayField = ['d' dayTokens{1}{1}];  % e.g., 'd20240316'
    
    % List subfolders in this day folder.
    subFolderInfo = dir(dayFolderPath);
    subFolderNames = {subFolderInfo([subFolderInfo.isdir]).name};
    % Remove '.' and '..'
    subFolderNames = setdiff(subFolderNames, {'.', '..'});
    
    for i = 1:length(subFolderNames)
        subFolderName = subFolderNames{i};
        subFolderPath = fullfile(dayFolderPath, subFolderName);
        
        % Parse subfolder name: expected format: subjectID_YYYYMMDDHHMMSS_task
        tokens = regexp(subFolderName, ['^' subjectID '_\d{14}_(\w+)$'], 'tokens');
        if isempty(tokens)
            warning('Subfolder "%s" does not match expected pattern.', subFolderName);
            continue;
        end
        taskType = lower(tokens{1}{1});  % e.g., 'training', 'decoding', etc.
        
        % Load the recording from this subfolder.
        taskSession = loadTaskData(subFolderPath);
        
        % Initialize the structure for this task type if needed.
        if ~isfield(tasks, taskType)
            tasks.(taskType) = struct();
        end
        
        % For this task type, add the session into the day grouping.
        if ~isfield(tasks.(taskType), dayField)
            tasks.(taskType).(dayField) = {taskSession};
        else
            tasks.(taskType).(dayField){end+1} = taskSession;
        end
    end
end

% Now, create the final data struct.
% For each task type, sort the day fields (they are named with a leading letter)
% and number them sequentially.
data = struct();
taskTypes = fieldnames(tasks);
for t = 1:length(taskTypes)
    taskType = taskTypes{t};  % e.g., 'decoding'
    dayFields = fieldnames(tasks.(taskType));
    dayFields = sort(dayFields); % sorting works because field names like 'd20240316'
    for i = 1:length(dayFields)
        combined = [];  % will hold concatenated data for this day and task
        sessions = tasks.(taskType).(dayFields{i});
        for j = 1:length(sessions)
            if isempty(combined)
                combined = sessions{j};
            else
                combined = concatenateSession(combined, sessions{j});
            end
        end
        % Create field name such as 'decoding1', 'training1', etc.
        fieldName = [taskType num2str(i)];
        data.(fieldName) = combined;
    end
end

end

%---------------------------------------------------------------------
function taskData = loadTaskData(folder)
% loadTaskData loads a recording session from a given folder.
%
% This example loads the first GDF file found in the folder.
files = dir(fullfile(folder, '*.gdf'));
if isempty(files)
    error('No GDF files found in folder %s', folder);
end
filePath = fullfile(files(1).folder, files(1).name);
[signal, header] = sload(filePath);
%
% Here we mimic your previous concatenation approach:
taskData.data = signal;
taskData.header = header;
% For example, you might record the end-of-file index as the number of samples.
taskData.eof = size(signal, 1);
end

%---------------------------------------------------------------------
function combined = concatenateSession(combined, newSession)
% concatenateSession combines two session structs that have fields:
%   data, header, and eof.
%
% This code assumes the following:
%  - Data matrices are concatenated along the first (time) dimension.
%  - Header.EVENT fields exist and need to be updated.
%
% If 'combined' is empty, simply return newSession.
if isempty(combined)
    combined = newSession;
    return;
end

% Offset: number of samples in the already combined data.
offset = size(combined.data, 1);

% Concatenate header EVENT fields if they exist.
if isfield(combined.header, 'EVENT') && isfield(newSession.header, 'EVENT')
    % Adjust new session event positions by the offset.
    newEventPos = newSession.header.EVENT.POS + offset;
    combined.header.EVENT.TYP = cat(1, combined.header.EVENT.TYP, newSession.header.EVENT.TYP);
    combined.header.EVENT.POS = cat(1, combined.header.EVENT.POS, newEventPos);
end

% Concatenate the data.
combined.data = cat(1, combined.data, newSession.data);

% Update the eof field: append the new end-of-data index.
combined.eof = cat(1, combined.eof, size(combined.data, 1));
end
