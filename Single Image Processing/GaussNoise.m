function [NewImage, Hnew] = GaussNoise(Image, Color, mean, var)
%UNTITLED5 Summary of this function goes here
%% 2) Gaussian Noise of user-specified parameters
ImageCopy = Image(:,:,Color);
NewImage = ImageCopy;
[W,L] = size(ImageCopy);
maxInt = 255;

for i = 1:W
    for j = 1:L
        noise = mean + randn(1)*sqrt(var);
        NewImage(i,j) = ImageCopy(i,j) + noise; 
    end
end

% Histogram Calculations
for i = 0:maxInt
    Hnew(i+1) = sum(sum(ImageCopy == i));
end

end

