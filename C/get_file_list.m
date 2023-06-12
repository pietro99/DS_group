function fileList = get_file_list(base_folder, folder)
    % Get a list of all subfolders within the input folder
    subfolders = dir(fullfile(base_folder, folder));
    subfolders = subfolders([subfolders.isdir]);
    subfolders = subfolders(3:end);  % Exclude "." and ".." directories
    
    fileList = [];
    % Loop through each subfolder
    for i = 1:numel(subfolders)
        subfolderPath = fullfile(base_folder,folder, subfolders(i).name);
        % Get a list of all image files within the subfolder
        imageFiles = dir(fullfile(subfolderPath, '*')); 
        imageFiles = imageFiles(3:end);
        fileList = [fileList; imageFiles];
    end
end