function [NewImage, Hnew] = MedFilter(Image, weight)
%ImgFilter: Process Image through a filter given weight inputs from user
% Median Filter
NewImage = Image; % copy original Image
[W, L] = size(NewImage); % width and length of image
inputs = []; % empty vector
Psum = 0; % sum of pixels
maxInt = 255; % max pixel intensity

w = size(weight, 1); % size of the weight matrix
n = (w + 1)/2; % indexing for the spatial neighborhood
c = n - 1; % move the position of the weights to expected MATLAB indexing
for u = n:W-n
    for v = n:L-n
        inputs = [];
        for i = -c:c
            for j = -c:c
                for val = 1:weight(i+n, j+n)   
                    inputs = [inputs, Image(u+i,v+j)];
                end
            end
        end
    NewImage(u,v) = median(inputs);       
    end
end

% Histogram Calculations
for i = 0:maxInt
    Hnew(i+1) = sum(sum(NewImage == i));
end

end

