%% Image Processing

%% 1) General framework for processing all images in a batch setting with supplying set-up initializing file
% Open image folder and add to search path
Folder = 'Cancerous cell smears';
addpath(Folder);
mainFolder = dir(Folder); % Stores all image names in struct format

name = 'Snake.JPG';
OImage = imread(mainFolder(3).name);
Image = OImage(:,:,1);

[X,Y] = size(Image);
subplot(2,3,1); imshow(Image);

[~, threshold] = edge(Image, 'sobel');
BWs = edge(Image, 'sobel', threshold*0.5);
subplot(2,3,2); imshow(BWs);

% Sobel Edge Detection
subplot(2,3,1); imshow(Image);
Gblur = [1 2 1; 2 4 2; 1 2 1];
Gscale = sum(sum(Gblur));
Hx = [-1 0 1; -2 0 2; -1 0 1];
Hy = [-1 -2 -1; 0 0 0; 1 2 1];
Scale = 1; 

Image = double(Image);
[ImY, HistY] = LinearFilter(Image, Hy, Scale); 
subplot(2,3,3); imshow(ImY);

[ImX, HistX] = LinearFilter(Image, Hx, Scale);
subplot(2,3,4); imshow(ImX);

Magnitude = sqrt(abs(ImX).^2 + abs(ImY).^2); 
ImComb = (Magnitude./max(max(Magnitude))).*255;
subplot(2,3,5); imshow(uint8(Magnitude));

% Convert to  binary and implement erosion and dilation operators
[W,L] = size(Image);
Magnitude = uint8(Magnitude);
subplot(2,2,1); imshow(Magnitude);
Imbin = zeros(W,L);
for i = 1:W
    for j = 1:L
        if Magnitude(i,j) >= 40
            Imbin(i,j) = 1;
        end
    end
end

subplot(2,2,2); imshow(Imbin);

% Dilate
B = [0 0; 1 1; -1 1; 1 -1; -1 -1]; % struct element
[ix, iy, v] = find(Imbin == 1);
newDil = false(W,L);
for n = 1:size(B,1)
    newMask = zeros(W,L);
    for i = 1: length(ix)
        newMask(ix(i) + B(n,1), iy(i) + B(n,2)) = Imbin(ix(i), iy(i));
    end
    newDil = newDil + newMask;
end

subplot(2,2,3); imshow(newDil);

% Erosion
newErd = zeros(W,L);
B = [0 1 0; 1 0 1; 0 1 0];

c=floor(size(B,1)/2);
d=floor(size(B,2)/2);

%Intialize a matrix with size of matrix A
for i=1:W-(2*c)
    for j=1:L-(2*d)     
        newMask = Imbin(i:i+(2*c),j:j+(2*d));      
        Ediff = min(min(newMask - B));
        if Ediff == 0
            newErd(i,j) = 1;
        end
    end
end
subplot(2,2,4); imshow(newErd);

%% Canny Edge Detection
% 1) Smooth image with Gaussian 
Fgauss = [1 2 1; 2 4 2; 1 2 1];
Fscale = sum(sum(Fgauss));

NImage = LinearFilter(Image, Fgauss, Fscale);

% 2) Compute Gradient magnitude and direction
Hx = [-1 0 1; -2 0 2; -1 0 1];
Hy = [-1 -2 -1; 0 0 0; 1 2 1];
Scale = 1; 

Image = double(Image);
[ImY, HistY] = LinearFilter(Image, Hy, Scale); 
subplot(2,3,3); imshow(ImY);

[ImX, HistX] = LinearFilter(Image, Hx, Scale);
subplot(2,3,4); imshow(ImX);

Magnitude = sqrt(abs(ImX).^2 + abs(ImY).^2); 
Direct = atan2(ImX, ImY);
Direct = Direct*180/pi;

subplot(2,3,5); imshow(uint8(Magnitude));

% 3) Thin edges by applying non-max suppression 

% 4) Detect edges by double thresholding

%% Image Segmentation
% 1) Histogram thresholding
Hist = Histogram(Image);
P = (Hist ./ (W*L))'; % probability of histogram intensities
Int = 0:255;
cdf = cumsum(P);

Thresh = 1;
maxV = 0;

for T = 2:255
    obj = Int(1:T);
    back = Int(T+1:end);
    Ub(T) = sum(back.*P(T+1:end)) / (1-cdf(T));
    Uo(T) = sum(obj.*P(1:T)) / cdf(T);
    VarB(T) = cdf(T)*(1-cdf(T))*(Uo(T) - Ub(T)).^2;
    
    if VarB(T) >= maxV
        maxV = VarB(T);
        Thresh = T;
    end
end

for i = 1:W
    for j = 1:L
        if Image(i,j) <= Thresh % Cell
            SegImg(i,j) = 1;
        elseif Image(i,j) > Thresh % Background
            SegImg(i,j) = 0;
        end        
    end
end

subplot(2,2,1); imshow(uint8(Image)); 

subplot(2,2,2); 
imshow(SegImg);

%% 2) K mean clustering
% Integrate color channel intensities and location  F(r,g,b,x,y)
Data = zeros(W*L, 1);
Red = OImage(:,:,1);
Data = double(Red(:));

x = [];
y = [];
for i = 1:W
    x = [x, ones(1,L)*i ];    
end
for i = 1:L
    y = [y, 1:W];
end

k = 5; % number of clusters
centers = 1:k;
Kmeans = zeros(1,k);
thresh = round((W*L) / k);
n = 1:thresh:(W*L);
n(k+1) = W*L;
flag = 0; 
t = 1;

for i = 1:k
    Kmeans(i) = mean(Data(n(i:i+1))); % mean of RBG for each center
end

KmeanNew = Kmeans;

while flag == 0
    Kmeans = KmeanNew;
    for i = 1:(W*L)
        for nn = 1:k
            Diff(nn) = sqrt((Kmeans(nn) - Data(i)).^2);
        end
        [val, idx(i)] = min(Diff);       
    end

    for nn = 1:k
        cluster = find(idx == nn);
        KmeanNew(nn) = mean(Data(cluster));
        Comp{nn} = cluster;
    end
    
    if sum(Kmeans - KmeanNew) == 0 || t > 5 % no change in mean values
        flag = 1;
    end
    
    t = t + 1;
    disp(['T = ', num2str(t)]);
    disp(['Kmean =', num2str(Kmeans)]);
end

% Making the new image based on the clusters found
NewImage = zeros(W,L);
Ints = linspace(0, 255, k+1);
Ints = round(Ints);

for i = 1:k
    idx = Comp{i};
    NewImage(idx) = Ints(i);
end

imshow(uint8(NewImage));