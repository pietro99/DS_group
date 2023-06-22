% Specify the input 
base_data_folder = 'data';
inputFolder = 'chest-rays';
outputFolder = inputFolder+"_processed_full_half";

process(base_data_folder,inputFolder, outputFolder)


function process(base_data_folder,inputFolder, outputFolder)
    % Create the data folder if it doesn't exist
    if ~isfolder(base_data_folder)
        mkdir(base_data_folder);
    end

    % Create the output folder if it doesn't exist
    if ~isfolder(fullfile(base_data_folder,outputFolder))
        mkdir(fullfile(base_data_folder,outputFolder));
    end
    
    % Get a list of all subfolders within the input folder
    subfolders = dir(fullfile(base_data_folder,inputFolder));
    subfolders = subfolders([subfolders.isdir]);
    subfolders = subfolders(3:end);  % Exclude "." and ".." directories
    
    % Loop through each subfolder
    for i = 1:numel(subfolders)
        subfolderPath = fullfile(base_data_folder, inputFolder, subfolders(i).name);
        % Get a list of all image files within the subfolder
        imageFiles = dir(fullfile(subfolderPath, '*')); 
        imageFiles = imageFiles(3:end);
        % Loop through each image file
        for j = 1:numel(imageFiles)
            imagePath = fullfile(subfolderPath, imageFiles(j).name);
            
            % Read the image
            image = imread(imagePath);
            
            % Perform image processing on the image (replace with your own processing code)
            processedImage = image_transform(image); % Example: resize the image
            
            % Generate the output file path and save the processed image
            [~, imageName, imageExt] = fileparts(imageFiles(j).name);
    
            new_subfolder = fullfile(base_data_folder, outputFolder, subfolders(i).name);
            if ~isfolder(new_subfolder)
                mkdir(new_subfolder);
            end
            outputFilePath = fullfile(new_subfolder, [imageName imageExt]);
            imwrite(processedImage, outputFilePath);
        end
    end
end

function processed_image = image_transform(image)
    image = im2gray(image);
    outputSize = [150, 150];
    processed_image = imresize(image, outputSize);
end

function fileList = get_file_list(folder)
fileList = get_file_list(folder);
end

