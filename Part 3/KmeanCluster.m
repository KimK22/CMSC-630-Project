function [NewImage, Comp] = KmeanCluster(Image, k)
% K-mean Cluster Function: Given initial image and K parameter, conduct K
% mean cluster analysis.
%%
centers = 1:k; % calculate center of clusters
[W,L] = size(Image);

% Set thresholds for the k parameter used (initial clusters and output
% intensities to separate clusters)
thresh = round((W*L) / k); 
indx = 1:thresh:(W*L);
indx(k+1) = W*L;
Ints = linspace(0, 255, k); % new output intensities
Ints = round(Ints);

t = 1;
flag = 0;

Data = double(Image(:));

Kmean = Ints;

KmeanNew = Kmean;

while flag == 0 
    % loop through k-means cluster algorithm until cluster means don't change
    Kmean = KmeanNew;
    for u = 1:(W*L)
        for v = 1:k
            Diff(v) = sqrt((Kmean(v) - Data(u)).^2); 
                % find euclidean distance between mean and data
        end
        [val, idx(u)] = min(Diff); 
            % assign to a cluster based on minimum difference      
    end

    for v = 1:k % calculate new cluster means
        cluster = find(idx == v);
        KmeanNew(v) = mean(Data(cluster));
        Comp{v} = cluster;
    end

    if sum(Kmean - KmeanNew) == 0 || t > 10 
        % if no change in mean values, stop looping
        flag = 1;
    end

    t = t + 1; % imposing restriction to prevent infinite looping
end

% Making the new image based on the clusters found
NewImage = zeros(W,L);

for v = 1:k
    idx = Comp{v};
    NewImage(idx) = Ints(v);
end

end

