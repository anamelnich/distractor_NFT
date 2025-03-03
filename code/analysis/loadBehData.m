function [behtraining, behvalidation, behdecoding] = loadBehData(path)
dirInfo = dir(path); %list all files in data folder
dirNames = {dirInfo.name}; %create a cell array containing names of all files and folders in data
dirNames = dirNames(~ismember(dirNames, {'.', '..','.DS_Store'}))
dirFolders = {dirInfo.folder}; %create a cell array of paths for each file and folder in data

trainingFile = contains(dirNames, 'training');
decodingFile = contains(dirNames, 'decoding');
validationFile = contains(dirNames, 'validation');
behtraining = []; behdecoding = []; behvalidation=[];
beh_var = {'trial', 'trial_type','response','tpos','dpos','dot','class'};
beh_var = {'trial', 'trial_type', 'response', 'tpos', 'dpos', 'dot', 'class'};
    for j = 1:length(beh_var)
        behtraining.(beh_var{j}) = [];
        behdecoding.(beh_var{j}) = [];
        behvalidation.(beh_var{j})=[];
    end
behtraining.runNumber = [];
behdecoding.runNumber = [];
behvalidation.runNumber = [];
behtraining.RT = [];
behdecoding.RT = [];
behvalidation.RT = [];

for i = 1:length(dirNames)
    exName = dirNames{i};
    exFolder = dirFolders{i};
    fileName = [exFolder '/' exName '/' exName];
    beh = readmatrix([fileName '.analysis.txt']);
    numColumns = size(beh,2);
    numTrials = size(beh, 1);
    runNumber = repelem(i, numTrials)';
    triggerData = readmatrix([fileName '.triggers.txt']);
    triggerData(triggerData(:, 2) == 6, :) = [];
    if trainingFile(i)
        for j = 1:numColumns
            behtraining.(beh_var{j}) = cat(1, behtraining.(beh_var{j}), beh(:, j));
        end
        behtraining.runNumber = [behtraining.runNumber; runNumber];
        trialStarts = triggerData(1:2:end, :); % Rows 1,3,5,...
        responses = triggerData(2:2:end, :);   % Rows 2,4,6,...
        RTs = responses(:, 3) - trialStarts(:, 3);
        behtraining.RT = [behtraining.RT; RTs];
    elseif validationFile(i)
        for j = 1:numColumns
            behvalidation.(beh_var{j}) = cat(1, behvalidation.(beh_var{j}), beh(:, j));
        end
         behvalidation.runNumber = [behvalidation.runNumber; runNumber];
        trialStarts = triggerData(1:2:end, :); % Rows 1,3,5,...
        responses = triggerData(2:2:end, :);   % Rows 2,4,6,...
        RTs = responses(:, 3) - trialStarts(:, 3);
        behvalidation.RT = [behvalidation.RT; RTs];
    elseif decodingFile(i)
        for j = 1:numColumns
            behdecoding.(beh_var{j}) = cat(1, behdecoding.(beh_var{j}), beh(:, j));
        end
         behdecoding.runNumber = [behdecoding.runNumber; runNumber];
        trialStarts = triggerData(1:2:end, :); % Rows 1,3,5,...
        responses = triggerData(2:2:end, :);   % Rows 2,4,6,...
        RTs = responses(:, 3) - trialStarts(:, 3);
        behdecoding.RT = [behdecoding.RT; RTs];
    end
end

