function [ImageEq, Heq] = HistEqual(Image, Color, H)
% HISTEQUAL: histogram equalization
% Goal: Compress based on peak values in the histogram 

ImageCopy = Image(:,:,Color);
ImageEq = ImageCopy;

[W, L] = size(ImageCopy);

maxInt = 255;
PixSum = 0;
CumSum = zeros(256,1);
FreqPix = zeros(256,1);
Heq = zeros(256,1);
PixTot = W*L;

for i = 1:maxInt + 1 % size(ProbVal)
    PixSum = PixSum + H(i);
    CumSum(i) = PixSum;
    FreqPix(i) = CumSum(i) / PixTot;
    Heq(i) = round(FreqPix(i)*maxInt);
end

for u = 1:W
    for v = 1:L
        ImageEq(u,v) = Heq(ImageCopy(u,v) + 1);
    end
end
