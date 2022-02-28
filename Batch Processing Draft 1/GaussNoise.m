function [NewImage, Hnew] = GaussNoise(Image, mean, var)
%UNTITLED5 Summary of this function goes here
%% 2) Gaussian Noise of user-specified parameters
NewImage = Image;
[W,L] = size(NewImage);
maxInt = 255;

for i = 1:W
    for j = 1:L
        noise = mean + randn(1)*sqrt(var);
        NewImage(i,j) = Image(i,j) + noise; 
    end
end

% Histogram Calculations
for i = 0:maxInt
    Hnew(i+1) = sum(sum(Image == i));
end

end

