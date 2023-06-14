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
L = strings(len, 1);
for i = 1:numel(fl)
    label = labels{i};
    for j = 1:numel(fl{i})
        image =imread(fullfile(fl{i}(j).folder , fl{i}(j).name));
        [h, w, d] = size(image);
        x = double(reshape(image, w*h*d, 1))/255;
        D(:, i) = x;
        L(j, 1) = label;
    end
end

num_of_rand = 10;
rand_pixels = randperm(size(D, 1),num_of_rand);
rand_images = randperm(size(D, 2),num_of_rand);

D_sampled_img =D(:,rand_images);
D_sampled_pixels =D(rand_pixels,:);
%PCA, Newdata = calc_mapping(10, D);

%%

Data_matrix = D;


%%

%constuct mean vector
D_means = mean(Data_matrix, 2);

%get average image
average_image =  uint8(reshape(D_means, h, w, d)*255);
%imshow(average_image)

%%

M = Data_matrix - D_means;

m = length(M);
d_length = length(Data_matrix);

image_reference =  uint8(reshape(Data_matrix(:,1), h, w, d)*255);
image_reference_diff =  uint8(reshape(M(:,1), h, w, d)*255);
%imshow(image_reference)
%imshow(image_reference_diff)

%%

%compute covariance matrix
C =  (1/m) * (M*M');

%G = (1/m) * (M'*M);


%%
%compute eigenvalues and eigenvectors
num_of_eigenvec = 10;
[Vec, D_val] = eigs(C, d);
eig_vals = diag(D_val);

%%
%compute basis vectors
U = Vec;

%U = (1 ./ sqrt(eig_vals))' .* (M * Vec);

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
%get eigenvectors weights for image
image_index = 3800;

new_dim = U' * Data_matrix;

weights = new_dim(:,image_index);

%show original image
imshow(reshape(D(:, image_index), h, w, d));

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
    new_image = new_image + weights_image(i)*U(:,i);
    new_image = rescale(new_image,0,1);
end
imshow(uint8(reshape(new_image, h, w, d)*255));


%%
function [PCA, Newdata] = calc_mapping(L, data)
    l=size(L,2);


    Index = 1 : size(data,2);

    FirstL = data(:,L);

    RMI = setdiff(Index,L);

    RestMat = data(:,RMI);

    Newdata = [FirstL RestMat];

    Newdata = Newdata - mean(Newdata,1);

    m = size(data,1);

    NyCov = (Newdata' * Newdata(:,1:l))./(m-1); %the natural way to compute it. 

    NyCovA = NyCov(1:l,1:l);

    NyCovB = NyCov((l+1):end,1:l);

    [eigvecA, eigvalA] =eig(NyCovA,'matrix');
    %computing the indices of the eigen values in decreasing order.
    [d,ind] = sort(diag(eigvalA),'descend');
    %sorting eigen values as per the computed indices.
    eigvalsorted = eigvalA(ind,ind);
    % similarly sorting the eigen vectors as per the computed indices.
    eigvecsorted = eigvecA(:,ind);

    U_A = eigvecsorted;

    U_B = NyCovB * U_A * inv(eigvalsorted); 

    PCA_Mat = [U_A; U_B];


    PCA = zeros(size(PCA_Mat));

    k = 1;   
    for a = L
        PCA(a,:)=PCA_Mat(k,:)/norm(PCA_Mat(k,:));
        k=k+1;
    end

    for a = RMI
        PCA(a,:)=PCA_Mat(k,:)/norm(PCA_Mat(k,:));
        k=k+1;
    end
end