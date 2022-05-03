%% General framework for processing all images in a batch setting with supplying set-up initializing file
% Open image folder and add to search path
Folder = 'Cancerous cell smears';
addpath(Folder);
mainFolder = dir(Folder); % Stores all image names in struct format
mainFolder(1:2) = [];
OImage = imread(mainFolder(58).name);

Orig = OImage(:,:,2);
Image = double(Orig);
[W,L] = size(Image);

% Set up Gradient Filters:
Gblur = [1 2 1; 2 4 2; 1 2 1];
Gscale = sum(sum(Gblur));
[NewImg2, Hn] = LinearFilter(Image, Gblur, Gscale);
% NewImg2 = Image - NewImg2;
[NewImg, clusters] = KmeanCluster(NewImg2, 5);

Imbin = zeros(W,L);
Nuclei = zeros(W,L);
for i = 1:W
    for j = 1:L
        if NewImg(i,j) < 190
            Imbin(i,j) = 1;
        end
    end
end

subplot(2,2,1); imshow(uint8(NewImg));
subplot(2,2,2); imshow(Imbin);

E = [0 1 0; 1 1 1; 0 1 0]; % erode image
D = [0 0; 1 1; -1 1; 1 -1; -1 -1]; % struct element

Eimg = ErodeImg(Imbin, E);
Eimg = ErodeImg(Eimg, E);
Eimg = ErodeImg(Eimg, E);
Dimg = DilateImg(Eimg, D);
subplot(2,2,3); imshow(Eimg);

Edge = Dimg - Eimg; 
Edge = DilateImg(Edge, D);
Edge = ErodeImg(Edge, E);

subplot(2,2,3); imshow(Edge);


ClassImg = zeros(W,L);
label = 1;
n = 4;
for i = n+1:W-n %n+1:W-n
    for j = n+1:L-n % n+1:L-n
         if Edge(i,j) == 1
            val = ClassImg(i-n:i+n, j-n:j+n);
            Oval = Edge(i-n:i+n, j-n:j+n);           
            if sum(sum(Oval)) > 0 && sum(sum(val)) == 0
              minV2 = label;
              label = label + 1;
            else
               ind = find(val(:) > 0);
               minV = val(ind);
               minV2 = min(minV);
            end
            
            ClassImg(i-n:i+n, j-n:j+n) = Oval.*minV2;
        end 
    end
end

cellnum = [];
for n = 1:label
    [c, d] = find(ClassImg == n);
    if length(c) > 500
        cellnum = [cellnum n];
    else
        ClassImg(c,d) = 0;
    end
end


Cell = randsample(cellnum, 1);
[c, d] = find(ClassImg == Cell); % Perimeter of cell
Perim = length(c);
CopyImg = zeros(W,L);
for i = 1:length(c)
    CopyImg(c(i), d(i)) = 1;    
end

% filling in masking area for cells
CellRad = 0;
CopyImg2 = CopyImg;
Jval = [];
for j = 1:L
    data = find(CopyImg(:,j) == 1);
    [ind2, ~] = max(data);
    [ind1, ~] = min(data);
    CopyImg2(ind1:ind2, j) = 1;
    
    if ((ind2 - ind1)/2) >= CellRad
        CellRad = (ind2 - ind1)/2;
        CellCenter = [j, CellRad + ind1];
    end
end

for i = 1:W
    data = find(CopyImg(i,:) == 1);
    [ind2, ~] = max(data);
    [ind1, ~] = min(data);
    CopyImg2(i, ind1:ind2) = 1;
end

CellPerim = sum(sum(CopyImg));
CellArea = sum(sum(CopyImg2));

subplot(2,2,4); imshow(CopyImg2); hold on; plot(CellCenter(1), CellCenter(2), 'r*');

% Masking to get average and max pixel intensity of the cell
maskImg = zeros(W,L); 
PixSum = 0;

for i = 1:W
    for j = 1:L
        if CopyImg2(i,j) == 1
            maskImg(i,j) = Image(i,j);
            PixSum = PixSum + Image(i,j); 
        end
    end
end

AvgInt = round(PixSum / CellArea);
maxInt = max(max(maskImg));