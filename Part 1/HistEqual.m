function [ImageEq, Heq] = HistEqual(Image, H)
% HISTEQUAL: histogram equalization
% Goal: Compress based on peak values in the histogram 

% Histogram Equalization
ImageEq = Image; % equalized image
[W, L] = size(ImageEq);
PixTot = W*L;

pSum = 0;
cdf = zeros(256, 1); % cdf
newInt = zeros(256, 1); % new pixel intensities
Heq = zeros(256,1);


for i = 1: length(H) 
    pSum = pSum + H(i); 
    cdf(i) = pSum / PixTot; % cdf calculation
    newInt(i) = floor(cdf(i)*maxInt); % new intensities
end

for u = 1:W
    for v = 1:L
        oldInt = Image(u,v); % find pixel value from original image
        ImageEq(u,v) = newInt(oldInt + 1); 
            % replace original value with new intensity
    end
end

for i = 0:maxInt
    Heq(i+1) = sum(sum(ImageEq == i)); % equalized histogram
end


end