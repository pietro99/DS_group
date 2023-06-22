%% Load and preprocess data

% Get a list of all subfolders

% folderPath = 'data/chest-rays_processed/';
% folderPath = 'data/Grapevine_Leaves_Image_Dataset_processed/';
folderPath = 'data/mnist_png_processed/';

subfolders = dir(folderPath);
subfolders = subfolders([subfolders.isdir]);
subfolders = subfolders(3:end);  % Exclude '.' and '..' directories

% extension = '*.jpeg';
extension = '*.png';

% Read the first image to get dimensions
firstSubfolder = fullfile(folderPath, subfolders(1).name);
images = dir(fullfile(firstSubfolder, extension));

firstImage = imread(fullfile(firstSubfolder, images(1).name));
[w, h, d] = size(firstImage);

% Initialize variables
numImages = 0;
for i = 1:numel(subfolders)
    if subfolders(i).isdir && ~strcmp(subfolders(i).name, '.') && ~strcmp(subfolders(i).name, '..')
        subfolderPath = fullfile(folderPath, subfolders(i).name);
        imageFiles = dir(fullfile(subfolderPath, extension)); 
        numImages = numImages + numel(imageFiles);
    end
end

% Initialize the imageVectors matrix
imageVectors = zeros(numImages, w*h*d);

% Loop over each subfolder and read images
index = 1;
images_nums = {}
for i = 1:numel(subfolders)
    subfolder = fullfile(folderPath, subfolders(i).name);
    images = dir(fullfile(subfolder, extension)); 
    images_nums = [images_nums, numel(images)];
    for j = 1:numel(images)
        imagePath = fullfile(subfolder, images(j).name);
        image = imread(imagePath);

        % Reshape and normalize the image, then save it in imageVectors
        imageVectors(index, :) = double(reshape(image, w * h * d, 1)) / 255;

        index = index + 1;
    end
end

%% Construct the squared distance matrix
% dissimilarity = zeros(size(imageVectors,1));

% for i=1:size(imageVectors,1)
%     for j =1:size(imageVectors,1)
%         dissimilarity(i,j) = pdist2(imageVectors(i,:),imageVectors(j,:),'squaredeuclidean');
%     end
% end
dissimilarity = squareform(pdist(imageVectors,'squaredeuclidean'));

%% Create the centering matrix and Gram matrix
n = size(dissimilarity,1);

identity = eye(n);

% one is qual to e * e'
one = ones(n);

% centering matrix
centering_matrix = identity - (1/n) * one;

Gram = -1/2 * centering_matrix * dissimilarity * centering_matrix;
%% Compute eigenvalues and vectors of 𝐺ram

% for any dimension ≤ n, choose 2 to visualize
dimension = 3;

[V,D] = eigs(Gram,dimension);

% MDS embeddings to a low-dimensional space.
X =  V * sqrt(D);

% project image onto eigenspace
px = imageVectors' * X * X' ;
px = px'

%% calculate the stress

dissimilarity_reduced = zeros(size(X,1));

for i=1:size(X,1)
    for j =1:size(X,1)
        dissimilarity_reduced(i,j) = pdist2(X(i,:),X(j,:),'squaredeuclidean');
    end
end

stress_matrix = (sqrt(dissimilarity_reduced) - sqrt(dissimilarity)) .* (sqrt(dissimilarity_reduced) - sqrt(dissimilarity))
% d_squre_sum = sum(dissimilarity_reduced(:))/2
stresses = sum(stress_matrix(:)) / 2;

%% test

%convert eigenvector to image and display the first image
new_images =uint8(reshape(rescale(px,0,1) , size(imageVectors,1),h,w,d)*255);
image1 = squeeze(new_images(1,:,:));

% Perform built-in MDS
[Y,eigvals] = cmdscale(sqrt(dissimilarity),dimension);

py = imageVectors' * Y * Y' ;
py = py'
new_images_builtin =uint8(reshape(rescale(py,0,1) , size(imageVectors,1),h,w,d)*255);
image1_builtin = squeeze(new_images_builtin(1,:,:));
% 
% dissimilarity_reduced_b = zeros(size(Y,1));
% 
% for i=1:size(Y,1)
%     for j =1:size(Y,1)
%         dissimilarity_reduced_b(i,j) = pdist2(Y(i,:),Y(j,:),'squaredeuclidean');
%     end
% end
% 
% stress_matrix_b = (sqrt(dissimilarity_reduced_b) - sqrt(dissimilarity)) .* (sqrt(dissimilarity_reduced_b) - sqrt(dissimilarity))
% d_squre_sum_b = sum(dissimilarity(:))/2
% stress_b = sqrt((sum(stress_matrix_b(:)) / 2)/d_squre_sum_b);


% Create a figure and divide it into a 1x3 grid of subplots
figure;
subplot(1, 3, 1);
imshow(image1);
title('MDS');

subplot(1, 3, 2);
imshow(image1_builtin );
title('built-in MDS');

subplot(1, 3, 3);
imshow(firstImage);
title('ORGINAL');

% figure;
% subplot(1, 2, 1);
% scatter(Y(:, 1), Y(:, 2));
% title('Bulit-in');
% 
% subplot(1, 2, 2);
% scatter(X(:, 1), X(:, 2));
% title('implementation');

%% scatter plot

% Assign class labels to each image
numImagesPerClass = cell2mat(images_nums);
classLabels = repelem(1:numel(numImagesPerClass), numImagesPerClass);

% Plot the scatter plot
figure;
hold on;

% Plot each class
for i = 1:numel(numImagesPerClass)
    classIndices = (classLabels == i);
    scatter(Y(classIndices, 1), Y(classIndices, 2), 'filled');
end

% Add labels and title
xlabel('x-axis');
ylabel('y-axis');
title('MDS - MNIST');

% Create a legend
legend(subfolders.name); % Update with appropriate class labels

% Hold off to finalize the plot
hold off;

%% draw line charts for stresses
% Dataset of dimensions
dimensions = [1, 2, 3, 4];

% Stresses for Chest X-ray
stresses_ChestXray = [3.1711e+07, 1.8478e+07, 1.2751e+07, 9.7383e+06];

% Stresses for MNIST
stresses_MNIST = [3.5329e+08, 2.4942e+08, 1.8830e+08, 1.4624e+08];

% Plotting the dataset for Chest X-ray
figure;
plot(dimensions, stresses_ChestXray, 'ro-', 'LineWidth', 2);
xlabel('Dimension');
ylabel('Stresses');
title('Stresses for Chest X-ray');
grid on;

% Plotting the dataset for MNIST
figure;
plot(dimensions, stresses_MNIST, 'bo-', 'LineWidth', 2);
xlabel('Dimension');
ylabel('Stresses');
title('Stresses for MNIST');
grid on;

