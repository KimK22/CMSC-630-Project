function [CellFeat, maskImg, maskImgC] = CellFeatExtract(Cbody, Coutline, Image, GImg)
%CellFeatExtract: Given cell isolation, find 5 different features
%% associated. 
[W,L] = size(Image);
CellPerim = sum(sum(Coutline));
CellArea = sum(sum(Cbody));

% Masking to get average and max pixel intensity of the cell
maskImg = zeros(W,L); 
maskImgC = 255*ones(W,L); 
PixSum = 0;

for i = 1:W
    for j = 1:L
        if Cbody(i,j) == 1
            maskImg(i,j) = Image(i,j);
            maskImgC(i,j) = GImg(i,j);
            
            PixInt = double(Image(i,j));
            
            PixSum = PixSum + PixInt; 
        end
    end
end

AvgInt = round(PixSum / CellArea);
MaxInt = max(max(maskImg));

[val, ind] = find(maskImgC < 80);
NucleiArea = length(ind); % Nuclei Area

CellFeat = [CellArea, CellPerim, AvgInt, MaxInt, NucleiArea];

end

