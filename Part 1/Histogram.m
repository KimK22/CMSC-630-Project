function [H] = Histogram(Image)
% HISTOGRAM: Calculate histogram for input image. 
% Histogram calculation for each original grayscale image
H = zeros(256,1);
maxInt = 255;

for i = 0:maxInt
    H(i+1) = sum(sum(Image == i));
end

end

