%specigy input folder
base_folder = 'data';
inputFolder = 'chest-rays';
outputFolder = inputFolder+"_processed";    

[fl, labels] = get_file_list(base_folder, outputFolder);

%%
%get size of images
example_image = imread(fullfile(fl{1}(1).folder , fl{1}(1).name));
[h, w, d] = size(example_image);
len = 0;
for i = 1:numel(fl)
    len = len + length(fl{i});
end
%construct data matrix 
D = zeros(w*h*d, len);
L = categorical(1,len);
counter = 0;
for i = 1:numel(fl)
    for j = 1:numel(fl{i})
        counter = counter +1;
        image =imread(fullfile(fl{i}(j).folder , fl{i}(j).name));
        [h, w, d] = size(image);
        x = double(reshape(image, w*h*d, 1))/255;
        D(:, counter) = x;
        L(1, counter) = labels{i};
    end
end
%%
num_of_rand = 1000;
rand_pixels = randperm(size(D, 1),num_of_rand);
rand_images = randperm(size(D, 2),num_of_rand);

D_sampled_img =D(:,rand_images);
D_sampled_pixels =D(rand_pixels,:);
%PCA, Newdata = calc_mapping(10, D);

%%

Data_matrix = D;
%%

U = nystrom(9000, D);


%%


%constuct mean vector
D_means = mean(Data_matrix, 2);

%get average image
average_image =  uint8(reshape(D_means, h, w, d)*255);
imshow(average_image)

%%
[h, w, d] = size(example_image);

M = Data_matrix - D_means;

m = length(M);
d_length = length(Data_matrix);

image_reference =  uint8(reshape(Data_matrix(:,1), h, w, d)*255);
image_reference_diff =  uint8(reshape(M(:,1), h, w, d)*255);
imshow(image_reference)
%imshow(image_reference_diff)

%%

%compute covariance matrix
C =  (1/m) * (M*M');

G = (1/m) * (M'*M);



%%
%compute eigenvalues and eigenvectors
num_of_eigenvec = 10;
[Vec, D_val] = eigs(C, num_of_eigenvec);
eig_vals = diag(D_val);

%%
%compute basis vectors
U = Vec;

%U = (1 ./ sqrt(eig_vals))' .* (M * Vec);

%%
%get eigenvectors weights for image
image_index = 311;

new_dim = U' * (M);

weights = new_dim(:,image_index);


%show original image
imshow(reshape(D(:, image_index), h, w, d));
scatter(new_dim(1,:), new_dim(2,:), 25, L, 'filled')
%scatter3(new_dim(1,:), new_dim(2,:),new_dim(3,:), 25, L, 'filled')

%%

%plot explained variance
explained_variance_tot = diag(D_val)/trace(C);
explained_variance_tot_sum = sum(explained_variance_tot);
explained_variance = diag(D_val)/trace(D_val);

figure;
plot(cumsum(explained_variance_tot))
ylim([0 1])
hold on
bar(explained_variance_tot)

%%
%plot eigenVectors as images
num1 = ceil(num_of_eigenvec/2);
num2 = 2;
clf
for i=1:num_of_eigenvec
    subplot(num1, num2, i);
    em = U(:, i);
    image = uint8(reshape(rescale(em, 0, 255), h, w, d));
    imshow(image);
end

%%

%plot projected eigenvectors of image
for i=1:5
    subplot(num1, num2, i);
    em = U(:, i);
    p1x= M(:, image_index)'*em*em;
    image = uint8(reshape(rescale(p1x, 0, 1), h, w, d)*255);
    imshow(image);
end

%%
%new_image = zeros(size(D_means));

%reconstruct image from eigenvectors combination
new_image = D_means;
%weights_image = [0,0,0,0,0]; if want to do it manually
for i = 1:length(U(1,:))
    new_image = new_image + weights(i)*U(:,i);
end
subplot(1, 2, 1);
imshow(uint8(reshape(new_image, h, w, d)*255));
subplot(1, 2, 2);
imshow(uint8(reshape(Data_matrix(:,image_index), h, w, d)*255));

%%
function [PCA] = nystrom(l,Data_matrix)
    Data_matrix = Data_matrix';
    [m, n] = size(Data_matrix);
    permutation_pixels = randperm(size(Data_matrix, 1));
    l_rand_pixels = permutation_pixels(1:l);

    rand_images_indeces = randperm(size(Data_matrix, 2));
    l_rand_images_indeces = rand_images_indeces(1:l);

    permutated_Data = Data_matrix(:,rand_images_indeces);
    
    %D_sampled_img =D(:,rand_images);
    %D_sampled_pixels =D(rand_pixels,:);

    estimated_C =  (1/m) * (permutated_Data' * permutated_Data(:,1:l));
    A = estimated_C(1:l, 1:l);
    B = estimated_C(l+1:end,1:l);

    [eig_vec, eig_val] = eigs(A, 10);

    UA = eig_vec;
    UB = B * UA * inv(eig_val); 
    PCA = [UA; UB];


end
