%specigy input folder
base_folder = 'data'
inputFolder = 'chest-rays';
outputFolder = inputFolder+"_processed";

fl = get_file_list(base_folder, outputFolder);

%%
example_image = imread(fullfile(fl(1).folder , fl(1).name));
[h, w, d] = size(example_image);
D = zeros(w*h*d, length(fl));
for i = 1:numel(fl)
    image =imread(fullfile(fl(i).folder , fl(i).name));
    [h, w, d] = size(image);
    x = double(reshape(image, w*h*d, 1))/255;
    D(:, i) = x;
end
%%

D_means = mean(D, 2);

average_image =  uint8(reshape(D_means, h, w, d)*255);
%imshow(average_image)
M = D - D_means;
m = length(M);
d_length = length(D);

image_reference =  uint8(reshape(D(:,1), h, w, d)*255);
image_reference_diff =  uint8(reshape(M(:,1), h, w, d)*255);
%imshow(image_reference)
%imshow(image_reference_diff)

%%

C =  (1/m) * (M*M');



%%
num_of_eigenvec = 10;
[Vec, D_val] = eigs(C, num_of_eigenvec);

%%

explained_variance_tot = diag(D_val)/trace(C);
explained_variance_tot_sum = sum(explained_variance_tot);
explained_variance = diag(D_val)/trace(D_val);
%%
figure;
plot(cumsum(explained_variance_tot))
ylim([0 1])

hold on

bar(explained_variance_tot)

%%
num1 = ceil(num_of_eigenvec/2);
num2 = 2;

%%
clf
for i=1:num_of_eigenvec
    subplot(num1, num2, i);
    em = Vec(:, i);
    image = uint8(reshape(rescale(em, 0, 255), h, w, d));
    imshow(image);
end
%%

image_index = 30900;

weights_image = Vec' * (M(:, image_index)- D_means);

imshow(reshape(D(:, image_index), h, w, d));

%%
for i=1:5
    subplot(num1, num2, i);
    em = Vec(:, i);
    p1x= M(:, image_index)'*em*em;
    image = uint8(reshape(rescale(p1x, 0, 1), h, w, d)*255);
    imshow(image);
end

%%

%image = uint8(reshape(p1x, h, w, d)*255);

%%

%figure, imshow(image)

%%

%u = M'*Vec;

%%
new_image = zeros(size(D_means));
new_image = D_means;
weights = [0,0,0,0,0];
for i = 1:length(Vec(1,:))
    new_image = new_image + weights_image(i)*Vec(:,i);
    new_image = rescale(new_image,0,1);
end
imshow(uint8(reshape(new_image, h, w, d)*255));