function behData = loadBehData(dataPath, subjectID)
% loadBehDataScalable loads behavioral data from multiple day folders.
%
% Assumptions:
%   - Day folders are named: subjectID_YYYYMMDD (e.g., 'e5_20240316')
%   - Within each day folder, subfolders are named: 
%         subjectID_YYYYMMDDHHMMSS_task
%       e.g., 'e5_20240316114601_training', 'e5_20240316120000_decoding',
%             'e5_20240316130000_stroop', or 'e5_20240316140000_EOGcalibration'
%
% For tasks 'training', 'decoding', or 'validation', sessions are loaded
% using the standard logic. For 'stroop', a different file (behoutput.txt)
% is loaded and the trial_type strings are converted. 'EOGcalibration'
% sessions are skipped.
%
% Output:
%   behData - a struct with fields like behData.training1, behData.decoding1, 
%             behData.stroop1, etc., where sessions from the same day 
%             (i.e., same day folder) are concatenated.

% List day folders for the subject (e.g., 'e5_20*')
dayFolders = dir(fullfile(dataPath, [subjectID '_20*']));
if isempty(dayFolders)
    error('No day folders found for subject %s in %s', subjectID, dataPath);
end

% Group sessions by task and day.
% We'll use a nested structure: tasks.(taskType).(dayField) = cell array of sessions.
tasks = struct();

for d = 1:length(dayFolders)
    dayFolderName = dayFolders(d).name;  % e.g., 'e5_20240316'
    dayFolderPath = fullfile(dayFolders(d).folder, dayFolderName);
    
    % Extract date from the day folder name.
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
        taskType = lower(tokens{1}{1});  % e.g., 'training', 'decoding', 'validation', 'stroop', or 'eogcalibration'
        
        % Skip EOGcalibration sessions.
        if strcmp(taskType, 'eogcalibration')
            continue;
        end
        
        % Load behavioral session using appropriate logic.
        behSession = loadBehSession(subFolderPath, taskType);
        disp(taskType)
        if isempty(behSession)
            continue;
        end
        
        % Group by task type and day.
        if ~isfield(tasks, taskType)
            tasks.(taskType) = struct();
        end
        if ~isfield(tasks.(taskType), dayField)
            tasks.(taskType).(dayField) = {behSession};
        else
            tasks.(taskType).(dayField){end+1} = behSession;
        end
    end
end

% Combine sessions for each task type and day.
behData = struct();
taskTypes = fieldnames(tasks);
for t = 1:length(taskTypes)
    taskType = taskTypes{t};  % e.g., 'training', 'decoding', or 'stroop'
    dayFields = fieldnames(tasks.(taskType));
    dayFields = sort(dayFields); % sorted by day, e.g., 'd20240316', 'd20240317',...
    for i = 1:length(dayFields)
        sessions = tasks.(taskType).(dayFields{i});
        combined = sessions{1};
        for j = 2:length(sessions)
            combined = concatenateBehSession(combined, sessions{j});
        end
        % Create field name such as 'training1', 'decoding1', 'stroop1', etc.
        fieldName = [taskType num2str(i)];
        behData.(fieldName) = combined;
    end
end

end

%---------------------------------------------------------------------
function behSession = loadBehSession(folder, taskType)
% loadBehSession loads behavioral data from a given session folder.
%
% If taskType is 'stroop', it loads the file ending in behoutput.txt, 
% skipping the first row, and converts the trial_type strings to integers:
%   neutral -> 0, congruent -> 1, incongruent -> 2.
%
% For other tasks ('training', 'decoding', 'validation'), it loads
% the .analysis.txt and .triggers.txt files.
%
% Input:
%   folder   - the path to the session folder
%   taskType - a string indicating the task type
%
% Output:
%   behSession - a struct containing the behavioral data

[~, sessionName] = fileparts(folder);

if strcmp(taskType, 'stroop')
    % Stroop logic: load behoutput.txt
    stroopFile = fullfile(folder, [sessionName, '.behoutput.txt']);
    if ~isfile(stroopFile)
        warning('Stroop file not found: %s', stroopFile);
        behSession = [];
        return;
    end
    % Set up import options to skip the first row (header) and use tab as delimiter
    opts = detectImportOptions(stroopFile, 'FileType', 'text', 'Delimiter', '\t');
    opts.DataLines = [2 Inf];  % start at second row
    T = readtable(stroopFile, opts);
    
    % Convert table to struct (each field is a column)
    behSession = table2struct(T, 'ToScalar', true);
    
    % Convert trial_type strings to integers
    % Assuming trial_type is stored as a cell array of strings.
    nTrials = height(T);
    new_trial_type = zeros(nTrials, 1);
    for i = 1:nTrials
        str = lower(T.Trial_Type{i});
        if strcmp(str, 'neutral')
            new_trial_type(i) = 0;
        elseif strcmp(str, 'congruent')
            new_trial_type(i) = 1;
        elseif strcmp(str, 'incongruent')
            new_trial_type(i) = 2;
        else
            new_trial_type(i) = NaN;
        end
    end
    behSession.trial_type = new_trial_type;
    % Optionally, you can rename or remove other fields if needed.
 
else
    % For training, decoding, or validation:
    analysisFile = fullfile(folder, [sessionName, '.analysis.txt']);
    triggersFile = fullfile(folder, [sessionName, '.triggers.txt']);
    
    if ~isfile(analysisFile)
        warning('Analysis file not found: %s', analysisFile);
        behSession = [];
        return;
    end
    if ~isfile(triggersFile)
        warning('Triggers file not found: %s', triggersFile);
        behSession = [];
        return;
    end
    
    beh = readmatrix(analysisFile);
    triggerData = readmatrix(triggersFile);
    % Remove rows with unwanted trigger values (e.g., 6)
    triggerData(triggerData(:,2) == 6, :) = [];

    
    % Define field names for standard behavioral variables.
    beh_vars = {'trial', 'trial_type', 'response', 'tpos', 'dpos', 'dot'};
    behSession = struct();
    for i = 1:length(beh_vars)
        behSession.(beh_vars{i}) = beh(:, i);
    end
    behSession.runNumber = [];  % Can be set later
    % Compute reaction times (RT) based on trigger events.
    if strcmp(taskType, 'decoding')  
        trialStarts = triggerData(triggerData(:,2)>50, :);
        %trialStarts = triggerData(1:3:end, :); % assume trial start events are in odd rows
        responses = triggerData(2:3:end, :);   % responses in even rows
    else
        trialStarts = triggerData(triggerData(:,2)>50, :); % assume trial start events are in odd rows
        responses = triggerData(2:2:end, :);   % responses in even rows
    end 
    RT = responses(:,3) - trialStarts(:,3);
    behSession.RT = RT;
end

end

%---------------------------------------------------------------------
function combined = concatenateBehSession(combined, newSession)
% concatenateBehSession combines two behavioral session structs.
%
% It concatenates the following fields:
%   - For standard tasks: 'trial', 'trial_type', 'response', 'tpos', 'dpos', 'dot', 'class'
%   - For stroop, the fields from the behoutput file (e.g., trial, trial_type, stimulus, ink_color, response, reaction_time)
%   - runNumber and RT fields are also concatenated.
%
% This function assumes the fields exist and are numeric (or arrays) and concatenates them vertically.
fields = fieldnames(newSession);
for i = 1:length(fields)
    field = fields{i};
    % If the field already exists in combined, concatenate vertically.
    if isfield(combined, field)
        combined.(field) = [combined.(field); newSession.(field)];
    else
        combined.(field) = newSession.(field);
    end
end
end
