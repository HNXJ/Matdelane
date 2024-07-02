function [idx, C] = jChannelClustering(x, numClusters, plotFlag)

%{
[idx, C] = jChannelClustering(x, numClusters, plotFlag)

>Input args 
x : Signal (N x M> numeric matrix, N-channels of length M)
numClusters : Cluster count (scalar, Hz)
plotFlag : Plotting status (boolean)

>Output args
idx : Assigned clusters (<1 x N> integer array)
C : Cluster centriods (<numClusters x M>  array) 

%}

    if ~exist("plotFlag", "var")

        plotFlag = 0;

    end

    disp("->K-means clustering with cosine metric on the x spectral response");

    if numClusters < 2

        numClusters = 2;
        warning("-->Number of the clusters cannot be less than 2");

    end
    
    imsx = zeros(size(x));

    for i = 1:size(x, 1)
    
        imsx(i, :) = abs(fft(x(i, :)));
    
    end

    [idx, C] = jCrossCluster(imsx, numClusters);

    if plotFlag

        figure;
        stem(idx);
        title("Cluster group ID");
        xlabel("Channels");
        ylabel("Groups");

    end

end
