function [H, pdf] = Histogram(Image, Color)
% HISTOGRAM: Calculate histogram for input image. 
% Histogram calculation for each original grayscale image
H = zeros(256,1);
maxInt = 255;
ImageCopy = Image(:,:,Color);
[W, L] = size(ImageCopy);

for i = 0:maxInt
    H(i+1) = sum(sum(ImageCopy == i));
    pdf(i+1) = H(i+1) / (W*L); % probability density (# pixels of intensity H / total # pixels)
end

end

