%% Load and preprocess data

% Specify the folder path
folderPath = 'data/chest-rays_processed/NORMAL';
% Get a list of all image files in the folder
imageFiles = dir(fullfile(folderPath, '*.jpeg'));

%get the size of images
example_image = imread(fullfile(folderPath, imageFiles(1).name));
[h, w, d] = size(example_image);

% Initialize a matrix to store the image vectors
imageVectors = zeros(length(imageFiles), h*w*d);

% Load and convert images
for i = 1:length(imageFiles)
    % Load the image
    imagePath = fullfile(folderPath, imageFiles(i).name);
    image = imread(imagePath);
    [h w d]=size(image);
    
    % Convert the image to a vector and store it in the cell array
    imageVectors(i,:) = double(reshape(image,w*h*d,1))/255;
end


%% Construct the similarity matrix
proximities = zeros(size(imageVectors,1));

for i=1:size(imageVectors,1)
    for j =1:size(imageVectors,1)
        proximities(i,j) = pdist2(imageVectors(i,:),imageVectors(j,:),'euclidean');
    end
end

%% Create the centering matrix and Gram matrix
n = size(proximities,1);

identity = eye(n);

% This is an equally sized matrix of 1s
one = ones(n);

% centering matrix
centering_matrix = identity - (1/n) * one;

Gram = -1/2 * centering_matrix * (proximities).*(proximities) * centering_matrix;
%% Compute eigenvalues and vectors of 𝐺ram

% for any dimension ≤ n, choose 2 to visualize
dimension = 2;

[V,D] = eigs(Gram,dimension);

%extract first eigenvector from matrix of eigenvectors
%ems=V(:,1);

% MDS embeddings to a low-dimensional space.
X =  V * sqrt(D);

%project image onto eigenspace
px = imageVectors' * X * X' ;

%convert eigenvector to image and display the first image
new_images =uint8(reshape(rescale(px',0,1) , length(imageFiles),h,w,d)*255);
imshow(squeeze(new_images(1,:,:)));
%imshow(new_images(1,:,:));

%% calculate the stress

proximities_reduced = zeros(size(new_images,1));

for i=1:size(new_images,1)
    for j =1:size(new_images,1)
        proximities_reduced(i,j) = pdist2(new_images(i,:),new_images(j,:),'euclidean');
    end
end

stress = (proximities_reduced - proximities)^2
