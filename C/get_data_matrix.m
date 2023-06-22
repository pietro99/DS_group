function [D,L, h, w, d] = get_data_matrix(fl, labels)
%get size of images
    example_image = imread(fullfile(fl{1}(1).folder , fl{1}(1).name));
    [h, w, d] = size(example_image);
    len = 0;
    for i = 1:numel(fl)
        len = len + length(fl{i});
    end
    %construct data matrix 
    D = zeros(w*h*d,len);
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