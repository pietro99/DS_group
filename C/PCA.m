%specigy input folder
base_folder = 'data';
inputFolder = 'chest-rays';
outputFolder = inputFolder+"_processed";

fl = get_file_list(base_folder, outputFolder);

%%
%get size of images
example_image = imread(fullfile(fl(1).folder , fl(1).name));
[h, w, d] = size(example_image);

%construct data matrix 
D = zeros(w*h*d, length(fl));
for i = 1:numel(fl)
    image =imread(fullfile(fl(i).folder , fl(i).name));
    [h, w, d] = size(image);
    x = double(reshape(image, w*h*d, 1))/255;
    D(:, i) = x;
end
%%

%constuct mean vector
D_means = mean(D, 2);

%get average image
average_image =  uint8(reshape(D_means, h, w, d)*255);
%imshow(average_image)

%%

M = D - D_means;
m = length(M);
d_length = length(D);

image_reference =  uint8(reshape(D(:,1), h, w, d)*255);
image_reference_diff =  uint8(reshape(M(:,1), h, w, d)*255);
%imshow(image_reference)
%imshow(image_reference_diff)

%%

%compute covariance matrix
C =  (1/m) * (M*M');

%%
%compute eigenvalues and eigenvectors
num_of_eigenvec = 10;
[Vec, D_val] = eigs(C, num_of_eigenvec);

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
    em = Vec(:, i);
    image = uint8(reshape(rescale(em, 0, 255), h, w, d));
    imshow(image);
end
%%
%get eigenvectors weights for image
image_index = 3900;

weights_image = Vec' * (M(:, image_index)- D_means);

%show original image
imshow(reshape(D(:, image_index), h, w, d));

%%

%plot projected eigenvectors of image
for i=1:5
    subplot(num1, num2, i);
    em = Vec(:, i);
    p1x= M(:, image_index)'*em*em;
    image = uint8(reshape(rescale(p1x, 0, 1), h, w, d)*255);
    imshow(image);
end

%%
%new_image = zeros(size(D_means));

%reconstruct image from eigenvectors combination
new_image = D_means;
%weights_image = [0,0,0,0,0]; if want to do it manually
for i = 1:length(Vec(1,:))
    new_image = new_image + weights_image(i)*Vec(:,i);
    new_image = rescale(new_image,0,1);
end
imshow(uint8(reshape(new_image, h, w, d)*255));