% Imports GDF file and creates Calibration, Training, and Decoding
% structure variables. 
% Calibration contains data, start index, end index, and header.
% Training contatins data, header (contains triggers), and end-of-file.
% Decoding contatins data, header (contains triggers), and end-of-file.
function [calibration, training, decoding] = loadData(path)
dirInfo = dir(path); %list all files in data folder
dirNames = {dirInfo.name}; %create a cell array containing names of all files and folders in data
dirNames = dirNames(~ismember(dirNames, {'.', '..','.DS_Store'}))
dirFolders = {dirInfo.folder}; %create a cell array of paths for each file and folder in data

calibrationFile = contains(dirNames, 'calibration'); %logical array where 1 indicates file name contains 'calibration'
trainingFile = contains(dirNames, 'training');
decodingFile = contains(dirNames, 'decoding');
calibration = []; training = []; decoding = []; 

for i = 1:length(dirNames)
    exName = dirNames{i};
    exFolder = dirFolders{i};
    fileName = [exFolder '/' exName '/' exName];
    [signal, header] = sload([fileName '.gdf']);
    
    if (calibrationFile(i))
        calibration.data = signal;
        calibration.start_index = find(signal(:, end) == 1);
        calibration.end_index = find(signal(:, end) == 2);
        calibration.header = header;
    elseif (trainingFile(i))
        if (isfield(training, 'data'))
            training.header.EVENT.TYP = cat(1, training.header.EVENT.TYP, header.EVENT.TYP);
            training.header.EVENT.POS = cat(1, training.header.EVENT.POS, header.EVENT.POS + size(training.data, 1));
            training.data = cat(1, training.data, signal);
            training.eof = cat(1, training.eof, size(training.data, 1));
        else
            training.data = signal;
            training.header = header;
            training.eof = size(signal, 1);
        end
    elseif (decodingFile(i))
        if (isfield(decoding, 'data'))
            decoding.header.EVENT.TYP = cat(1, decoding.header.EVENT.TYP, header.EVENT.TYP);
            decoding.header.EVENT.POS = cat(1, decoding.header.EVENT.POS, header.EVENT.POS + size(decoding.data, 1));
            decoding.data = cat(1, decoding.data, signal);
            decoding.eof = cat(1, decoding.eof, size(decoding.data, 1));
        else
            decoding.data = signal;
            decoding.header = header;
            decoding.eof = size(signal, 1);
        end
    end
end
delete sopen.mat
end
