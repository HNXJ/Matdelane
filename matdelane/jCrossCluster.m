function [idx, C] = jCrossCluster(x, numClusters)

%{

[idx, C] = jCrossCluster(x, numClusters)

>Input args 
x : Signal (N x M> numeric matrix, N-channels of length M)
numClusters : Cluster count (scalar, Hz)

>Output args
idx : Assigned clusters (<1 x N> integer array)
C : Cluster centriods (<numClusters x M>  array) 

%}

    [N, T] = size(x);

    if T < 2

        error("Signal length T must be greater than 1.");

    end

    similarityMatrix = zeros(N, N);

    parfor i = 1:N

        for j = 1:N

            corrValue = max(xcorr(x(i, :), x(j, :), 'coeff'));
            similarityMatrix(i, j) = corrValue;

        end

    end

    similarityMatrix(isnan(similarityMatrix)) = 0;
    distanceMatrix = 1 - similarityMatrix;

    opts = statset('UseParallel', 1);
    [idx, C] = kmeans(distanceMatrix, numClusters, 'Distance', 'cosine', 'Options', opts);

end
