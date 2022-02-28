function [NewImage, Hnew, MSQE] = ImageQuant(Image, Color, Level)
%ImageQuant: 
% Goal: Reduce number of bins in histogram / reduce values a pixel can be
% in the image.
%%
ImageCopy = Image(:,:,Color);
NewImage = ImageCopy;
[W, L] = size(ImageCopy);

maxInt = 255;
spaces = floor(maxInt / Level);
Thresholds = 0:spaces:maxInt;

PixInt = zeros(1,Level);

for nn = 2:Level
    PixVals = [Thresholds(nn-1):Thresholds(nn)];
    TotSum = sum(PixVals);
    total = length(PixVals);
    PixInt(nn) = TotSum/total;
end

% Replace pixel intensity values in the new image copy
for u = 1:W
    for v = 1:L
        for nn = 1:Level-1
            if ImageCopy(u,v) >= Thresholds(nn) && ImageCopy(u,v) ...
                    < Thresholds(nn+1)
                NewImage(u,v) = PixInt(nn);
            end
        end
    end
end

% Histogram calculation for quantized image
for i = 0:maxInt
    Hnew(i+1) = sum(sum(NewImage == i));
end

% MSQE for image quantization levels
MSQE = sum((ImageCopy(:) - NewImage(:)).^2) / (W*L);

end

