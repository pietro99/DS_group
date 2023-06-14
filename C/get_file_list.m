function [fileList, labels] = get_file_list(base_folder, folder)
    % Get a list of all subfolders within the input folder
    subfolders = dir(fullfile(base_folder, folder));
    subfolders = subfolders([subfolders.isdir]);
    subfolders = subfolders(3:end);  % Exclude "." and ".." directories
    
    fileList = cell(numel(subfolders), 1);
    labels = cell(numel(subfolders), 1);
    labeledData = table();
    % Loop through each subfolder
    for i = 1:numel(subfolders)
        subfolderPath = fullfile(base_folder,folder, subfolders(i).name);
        % Get a list of all image files within the subfolder
        imageFiles = dir(fullfile(subfolderPath, '*')); 
         
        imageFiles = imageFiles(3:end);
        fileList{i} = imageFiles;
        labels{i} = subfolders(i).name;
        
        %fileList = [fileList; imageFiles];
    end
end