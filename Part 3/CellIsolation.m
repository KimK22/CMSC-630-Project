function [CellArea, CellPerim] = CellIsolation(Edge)
%CELL ISOLATION: If the total cell area in the image is less than 1/2 of
%the whole image, identify one cell in the image to feature extract from. 

[W,L] = size(Edge); 
NewImg = zeros(W,L);
label = 1; % labelling the number of unique contours in the image
n = 4; % window size

for i = n+1:W-n
    for j = n+1:L-n
         if Edge(i,j) == 1
            val = NewImg(i-n:i+n, j-n:j+n); % see if a label was put in
            Oval = Edge(i-n:i+n, j-n:j+n); % see if an edge is present          
            if sum(sum(Oval)) > 0 && sum(sum(val)) == 0
              minV2 = label;
              label = label + 1; % set up new object identification
            else
               ind = find(val(:) > 0);
               minV = val(ind);
               minV2 = min(minV);
               if max(minV) - minV2 > 0
                   ind = find(NewImg == max(minV));
                   NewImg(ind) = minV2;
               end
            end

            NewImg(i-n:i+n, j-n:j+n) = Oval.*minV2; % replace values with new labels
        end 
    end
end


cellnum = []; % find number of actual cells in the image
for n = 1:label
    [c, d] = find(NewImg == n);
    if length(c) > 600
        cellnum = [cellnum n];
    else
        NewImg(c,d) = 0;
    end
end
 
% Chose one cell to save and conduct feature extraction on
Cell = randsample(cellnum, 1);
[c, d] = find(NewImg == Cell); % Perimeter of cell

CellPerim = zeros(W,L);
for i = 1:length(c)
    CellPerim(c(i), d(i)) = 1;    
end

% filling in masking area for cells
CellRad = 0;
CellArea = CellPerim;
Jval = [];
for j = 1:L
    data = find(CellPerim(:,j) == 1);
    [ind2, ~] = max(data);
    [ind1, ~] = min(data);
    CellArea(ind1:ind2, j) = 1;
end

for i = 1:W
    data = find(CellPerim(i,:) == 1);
    [ind2, ~] = max(data);
    [ind1, ~] = min(data);
    CellArea(i, ind1:ind2) = 1;
end

end
