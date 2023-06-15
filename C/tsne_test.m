load('mnist.mat');

rescaled = real(training.images(:,:,1:500)*255);
[h, w, d] = size(rescaled);
X = double(transpose(reshape(rescaled, w*h, d)));
labels = training.labels(1:500);

Y = run_tsne(X);
subplot(2,2,1)
gscatter(Y(:,1),Y(:,2), labels)
title('TSNE-MNIST')

Y_mds = run_mds(X);
subplot(2,2,3)
gscatter(Y_mds(:,1),Y_mds(:,2), labels)
title('MDS-MNIST')

[X_xray, xray_labels] = getXrayData();
X_xray = transpose(X_xray);
xray_labels = transpose(xray_labels);

X_slice = X_xray(1:500,:);
xray_labels = xray_labels(1:500,:);

Y_xray = run_tsne(X_slice);
subplot(2,2,2)
gscatter(Y_xray(:,1),Y_xray(:,2), xray_labels)
title('TSNE-XRAY')

Y_xray_mds = run_mds(X_slice);
subplot(2,2,4)
gscatter(Y_xray_mds(:,1),Y_xray_mds(:,2), xray_labels)
title('MDS-XRAY')

function Y = run_tsne(X)
    rng('default')
    Y = tsne(X,'Algorithm','exact','Distance','cosine');
end

function Y_mds = run_mds(X)
    rng('default')
    D = pdist(X);
    Y_mds = mdscale(D,2,'Start', 'random');
end

function [D,L] = getXrayData()
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
end