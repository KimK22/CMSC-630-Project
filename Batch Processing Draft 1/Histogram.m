function [H, pdf] = Histogram(Image)
% HISTOGRAM: Calculate histogram for input image. 
% Histogram calculation for each original grayscale image
H = zeros(256,1);
maxInt = 255;
[W, L] = size(Image);

for i = 0:maxInt
    H(i+1) = sum(sum(Image == i));
    pdf(i+1) = H(i+1) / (W*L); % probability density (# pixels of intensity H / total # pixels)
end

end

