function file_paths = select_files()
    % Opens a file dialog for the user to select multiple text files.
    %
    % Returns:
    %   file_paths (cell array of strings): Paths to the selected text files.
    
    [files, path] = uigetfile('*.txt', 'Select Text Files', '../', 'MultiSelect', 'on');
    
    if isequal(files, 0)
        % User clicked Cancel
        file_paths = {};
        disp('No files selected.');
    else
        if ischar(files)
            % If only one file is selected, uigetfile returns a character array
            files = {files};
        end
        % Concatenate path and filenames
        file_paths = fullfile(path, files);
    end
end