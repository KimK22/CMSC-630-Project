function [NewImage] = DilateImg(Image, D)
% Dilate: Dilate image given the structuring element provided
[W,L] = size(Image);
ImgCopy = zeros(W,L); 
ImgCopy(2:end-1, 2:end-1) = Image(2:end-1, 2:end-1);
[ix, iy, v] = find(ImgCopy == 1);
NewImage = zeros(W,L);
for n = 1:size(D,1)
    newMask = zeros(W,L);
  for c = 1:length(ix)
    newMask(ix(c) + D(n,1), iy(c) + D(n,2)) = ImgCopy(ix(c), iy(c));
  end
    NewImage = NewImage + newMask;
end

for i = 1:W
    for j = 1:L
        if NewImage(i,j) > 1
            NewImage(i,j) = 1;
        end
    end
end

end

