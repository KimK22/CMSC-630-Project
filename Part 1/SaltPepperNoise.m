function [NewImage, Hnew] = SaltPepperNoise(Image, Dnoise)
%SaltPepperNoise: Add noise to the image
NewImage = Image;
[W,L] = size(NewImage);
maxInt = 255;

Noise = rand(W, L);

for i = 1:W
    for j = 1:L
        if Noise(i,j) > 0 && Noise(i,j) < Dnoise / 2
            NewImage(i,j) = 0;
        elseif Noise(i,j) > Dnoise/2 && Noise(i,j) < Dnoise
            NewImage(i,j) = maxInt;
        end         
    end
end

% Histogram Calculation
for i = 0:maxInt
    Hnew(i+1) = sum(sum(NewImage == i));
end


end

