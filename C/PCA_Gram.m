%specigy input folder
base_folder = 'data';
inputFolder = 'MINIST';
outputFolder = inputFolder+"_processed_150";    

[fl, labels] = get_file_list(base_folder, inputFolder);
[D, L, h, w, d] = get_data_matrix(fl, labels);

%select subset
num_of_rand = 10000;
rand_images = randperm(size(D, 2),num_of_rand);

D_sampled_img =D(:,rand_images);
L = L(:, rand_images);


Data_matrix = D_sampled_img;




example_image = uint8(reshape(Data_matrix(:, 1), h, w, d));
[h, w, d] = size(example_image);

%constuct mean vector
D_means = mean(Data_matrix, 2);

M = Data_matrix - D_means;

m = length(M);
d_length = length(Data_matrix);

image_reference =  uint8(reshape(Data_matrix(:, 1), h, w, d)*255);
image_reference_diff =  uint8(reshape(M(:,1), h, w, d)*255);
%imshow(image_reference)
%imshow(image_reference_diff)


%compute gram matrix

G = (1/(m-1)) * (M'*M);



%compute eigenvalues and eigenvectors
num_of_eigenvec = 10;
[Vec, D_val] = eigs(G,num_of_eigenvec );
eig_vals = diag(D_val);

%compute basis vectors
U = (1 ./ sqrt(eig_vals))' .* (M * Vec);

%get eigenvectors weights for image
image_index = 311;

new_dim = U'*M;


%%
weights = new_dim(:,image_index);


%show original image
%imshow(reshape(Data_matrix(:, image_index), h, w, d));
scatter(new_dim(1,:), new_dim(2,:), 25, L, 'filled')
%scatter3(new_dim(1,:), new_dim(2,:),new_dim(3,:), 25, L, 'filled')

%%

%plot explained variance
explained_variance_tot = diag(D_val);
total_var = trace(C);
explained_variance_tot_sum = sum(explained_variance_tot);
explained_variance = diag(D_val);

figure;


plot(cumsum(explained_variance))
hold on

bar(explained_variance_tot)
ylim([0,total_var ])
xlabel('Principal Component')
ylabel('Variance')
yyaxis right

ylabel('Variance Precentage')
ax = gca;
ax.YAxis(2).Color = 'black';
title("Grapevine leaves")

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


projected_images = rescale(((new_dim' * U')+D_means')', 0,1);
distances = sqrt(sum((Data_matrix - projected_images).^2, 1));
[max_val, max_index] = max(distances);
[min_val, min_index] = min(distances);
save("distance-"+inputFolder+"-Gram.mat","distances")

subplot(1, 2, 1);

image_index = 300;

imshow(uint8(reshape(rescale(projected_images(:, max_index), 0, 1), h, w, d)*255));

subplot(1, 2, 2);
imshow(uint8(reshape(Data_matrix(:,max_index), h, w, d)*255));


%%
%%get image with manual weights
new_image = D_means;
weights = [0,0,0,0,0, 0, 0, 0, 0, 0]; 
for i = 1:length(U(1,:))
    new_image = new_image + weights(i)*U(:,i);
end
subplot(1, 2, 1);
imshow(uint8(reshape(new_image, h, w, d)*255));
subplot(1, 2, 2);
imshow(uint8(reshape(Data_matrix(:,image_index), h, w, d)*255));
%%
load("distance-chest-rays.mat", "distances");
d1 = distances;
load("distance-chest-rays-Gram.mat", "distances");
d2 = distances;
load("distance-chest-raysNystrom.mat", "distances");
d3 = distances;
load("distance-MINIST.mat", "distances");
d4 = distances;
load("distance-MINIST-Gram.mat", "distances");
d5 = distances;
load("distance-MINISTNystrom.mat", "distances");
d6 = distances;
load("distance-Grapevine_Leaves_Image_Dataset.mat", "distances");
d7 = distances;
load("distance-Grapevine_Leaves_Image_Dataset-Gram.mat", "distances");
d8 = distances;
load("distance-Grapevine_Leaves_Image_DatasetNystrom.mat", "distances");
d9 = distances;

dist= [d1, d2, d3, d4,d5,d6, d7, d8, d9];
values = {'chest-rays','chest-rays Snapshot','chest-rays Nystrom', 'MINIST','MINIST Snapshot','MINIST Nystrom', 'Grapevine Leaves',  'Grapevine Leaves w/Gram','Grapevine Leaves Nystrom'};
lengths = [size(d1, 2), size(d2, 2), size(d3,2),size(d4, 2), size(d5, 2), size(d6,2),size(d7, 2), size(d8, 2), size(d9,2)];
% Create an empty color matrix
grp = categorical(repelem(values, lengths));

% Create categorical array
%grp = [zeros(size(d1)),ones(size(d2)), ones(size(d3))+1];
boxplot(dist', grp, color = boxplotColors);
ylabel("distance")


%%
function [PCA] = nystrom(l,Data_matrix)
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

