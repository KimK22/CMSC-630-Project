function [NewImage, Hnew, MSQE] = ImageQuant(Image, Level)
%ImageQuant: 
% Goal: Reduce number of bins in histogram / reduce values a pixel can be
% in the image.
%%
NewImage = Image;
[W, L] = size(Image);

maxInt = 255;
spaces = floor(maxInt / Level);
Thresh = 0:spaces:maxInt;

PixInt = zeros(1,Level);
Hnew   = zeros(1,256);

for nn = 2:Level
    PixVals = Thresh(nn-1):Thresh(nn);
    TotSum = sum(PixVals);
    total = length(PixVals);
    PixInt(nn) = TotSum/total;
end

% Replace pixel intensity values in the new image copy
for u = 1:W
    for v = 1:L
        for nn = 1:Level-1
            if Image(u,v) >= Thresh(nn) && Image(u,v) < Thresh(nn+1)
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
MSQE = sum((Image(:) - NewImage(:)).^2) / (W*L);

end

