function [NewImage, Hnew] = LinearFilter(Image, weight, scale)
%ImgFilter: Process Image through a filter given weight inputs from user
[W, L] = size(Image); % width and length of image
NewImage = zeros(W,L); % copy original Image
Psum = 0; % sum of pixels
maxInt = 255; % max pixel intensity

w = size(weight,1);
n = (w + 1)/2; % indexing for the spatial neighborhood
c = n - 1; % move the position of the weights to expected MATLAB indexing

if scale == 0 % if this is a difference vector
    scale = 1; % set scale to 1 so it doesn't divide by 0
end

% Linear Filter
for u = n:W-n
    for v = n:L-n
        Psum = 0; % sum of average
        for i = -c:c
            for j = -c:c
              Psum = Psum + ((Image(u+i, v+j)*weight(i+n, j+n))/scale);
            end
        end
        NewImage(u,v) = round(Psum);
    end
end

%NewImage = uint8(NewImage);

% Histogram Calculations
for i = 0:maxInt
    Hnew(i+1) = sum(sum(NewImage == i));
end

end
