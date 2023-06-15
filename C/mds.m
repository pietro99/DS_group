%% Load and preprocess data

% Specify the folder path
folderPath = 'data/chest-rays_processed/NORMAL/';
% folderPath = 'data/Pistachio_Image_Dataset_processed/Kirmizi_Pistachio/';

% Get a list of all image files in the folder
imageFiles = dir(fullfile(folderPath, '*.jpeg'));
% imageFiles = dir(fullfile(folderPath, '*.jpg'));

%get the size of images from the first image
example_image = imread(fullfile(folderPath, imageFiles(1).name));
[h, w, d] = size(example_image);

% Initialize a matrix to store the images, each row is an image vector
imageVectors = zeros(length(imageFiles), h*w*d);

% Load, convert and save the images
for i = 1:length(imageFiles)
    imagePath = fullfile(folderPath, imageFiles(i).name);
    image = imread(imagePath);
    [h w d]=size(image);
    % Convert the image to a vector and store it in the image matrix
    imageVectors(i,:) = double(reshape(image,w*h*d,1))/255;
end

%% Construct the similarity matrix, each element is pair-wise squared distance
similarity = zeros(size(imageVectors,1));

for i=1:size(imageVectors,1)
    for j =1:size(imageVectors,1)
        similarity(i,j) = pdist2(imageVectors(i,:),imageVectors(j,:),'squaredeuclidean');
    end
end

%% Create the centering matrix and Gram matrix
n = size(similarity,1);

identity = eye(n);

% one is qual to e * e'
one = ones(n);

% centering matrix
centering_matrix = identity - (1/n) * one;

Gram = -1/2 * centering_matrix * similarity * centering_matrix;
%% Compute eigenvalues and vectors of 𝐺ram

% for any dimension ≤ n, choose 2 to visualize
dimension = 4;

[V,D] = eigs(Gram,dimension);

% MDS embeddings to a low-dimensional space.
X =  V * sqrt(D);

% project image onto eigenspace
px = imageVectors' * X * X' ;
px = px'

%% calculate the stress

similarity_reduced = zeros(size(X,1));

for i=1:size(X,1)
    for j =1:size(X,1)
        similarity_reduced(i,j) = pdist2(X(i,:),X(j,:),'squaredeuclidean');
    end
end

stress_matrix = (sqrt(similarity_reduced) - sqrt(similarity)) .* (sqrt(similarity_reduced) - sqrt(similarity))
stress = sum(stress_matrix(:)) / 2;

%% test

%convert eigenvector to image and display the first image
new_images =uint8(reshape(rescale(px,0,1) , length(imageFiles),h,w,d)*255);
image1 = squeeze(new_images(1,:,:));

% Perform built-in MDS with 2 dimensions
[Y,stress_builtin] = mdscale(similarity, dimension,'Start', 'random');

py = imageVectors' * Y * Y' ;
py = py'
new_images_builtin =uint8(reshape(rescale(py,0,1) , length(imageFiles),h,w,d)*255);
image1_builtin = squeeze(new_images_builtin(1,:,:));

similarity_reduced_b = zeros(size(Y,1));

for i=1:size(Y,1)
    for j =1:size(Y,1)
        similarity_reduced_b(i,j) = pdist2(Y(i,:),Y(j,:),'squaredeuclidean');
    end
end

stress_matrix_b = (sqrt(similarity_reduced_b) - sqrt(similarity)) .* (sqrt(similarity_reduced_b) - sqrt(similarity))
stress_b = sum(stress_matrix_b(:)) / 2;


% Create a figure and divide it into a 1x3 grid of subplots
figure;
subplot(1, 3, 1);
imshow(image1);
title('MDS');

subplot(1, 3, 2);
imshow(image1_builtin );
title('built-in MDS');

subplot(1, 3, 3);
imshow(example_image);
title('ORGINAL');

% figure;
% subplot(1, 2, 1);
% scatter(Y(:, 1), Y(:, 2));
% title('Bulit-in');
% 
% subplot(1, 2, 2);
% scatter(X(:, 1), X(:, 2));
% title('implementation');
