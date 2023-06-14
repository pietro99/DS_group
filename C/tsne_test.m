load('mnist.mat');

function X = run_tsne()
    rescaled = real(training.images(:,:,1:500)*255);
    [h, w, d] = size(rescaled);
    X = double(transpose(reshape(rescaled, w*h, d)));
    Y = tsne(X,'Algorithm','exact','Distance','cosine');
    gscatter(Y(:,1),Y(:,2), training.labels(1:500))
    title('MNIST')
end

function X_mds = run_mds(D)
    X_mds = mdscale(D,2);
end