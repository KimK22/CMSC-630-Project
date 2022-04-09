% Name: Kristin Kim
% Course: CMSC 630
% Outline: Project Part 2

% Goal: 
%     Call 500 real-life cancerous smear images and edit them using
%     different functionalities. This will be the main setup file that has 
%     the folder/ images imported and runs the functions defined. 

% Requirements: 
%     1) Batch setting image processing given folder location
%     2) Implement 1 selected edge detection algorithm
%     3) Implement Erosion and Dilation Operators
%     4) Implement 2 segmentation Techniques
%           A) Histogram Thresholding: Divide image into cells and
%           background
%           B) Clustering: Look at effects of different k parameters on
%           segmentation
%     5) Display Performance Measurements: 
%             Processing time (per procedure and average per image)

% Note: Super045.BMP = corrupted file, so I removed it

%% Cleanup Section
clc; 
close all; 
clear all;

%% 1) Batch Image Initialization
Folder = 'Cancerous cell smears'; % Location of images
addpath(Folder);
mainFolder = dir(Folder); % Stores all image names in struct format
mainFolder(1:2) = [];
% Separate different classes into variables
for j = 1:499 % length(mainFolder)
    if mainFolder(j).name(1:3) == 'cyl'
        cyl{j} = imread(mainFolder(j).name);
    elseif mainFolder(j).name(1:3) == 'int'
        inter{j-50} = imread(mainFolder(j).name);
    elseif mainFolder(j).name(1:3) == 'let'
        let{j-100} = imread(mainFolder(j).name);
    elseif mainFolder(j).name(1:3) == 'mod'
        mod{j - 200} = imread(mainFolder(j).name);
    elseif mainFolder(j).name(1:3) == 'par'
        para{j - 300} = imread(mainFolder(j).name);
    elseif mainFolder(j).name(1:3) == 'sup'
        super{j - 350} = imread(mainFolder(j).name);
    elseif mainFolder(j).name(1:3) == 'sva'
        svar{j - 399} = imread(mainFolder(j).name);
    end
end

%% Setup Variables
%  Integrate categories into same cell
CellType  = {cyl, inter, let, mod, para, super, svar};
CellName  = {'cyl', 'inter', 'let', 'mod', 'para', 'super', 'svar'};

%  Choosing Color Channel
Color     = 1;      % Red = 1; Green = 2; Blue = 3

% Size of Image
[W,L] = size(cyl{1}(:,:,1));
%% 2) Create Class structure with all images, color conversion
Class = struct();
numClass = length(CellType);
Htime = zeros(1,7); 
avgTime = zeros(1,7);
Hsum = zeros(256,1);
Havg = zeros(256,7);

for i = 1:numClass
    tic
    Class(i).name  = CellName{i};
    CellColor = {};
    H = {};
    Hsum = zeros(256,1);
    numImg = length(CellType{i});
    for j = 1:numImg
        Ind = cell2mat(CellType{i}(j));
        CellColor{j} = Ind(:,:,Color);
        H{j} = Histogram(CellColor{j});
        Hsum = Hsum + H{j}; 
    end
    
    Class(i).image = CellColor;
    Class(i).hist  = H;
    
    Havg(:,i) = round(Hsum ./ numImg);
    Htime(i) = toc;
    avgTime(i) = Htime(i) / numImg;
    disp(['Total Processing Time (s): ', num2str(Htime(i))]);
    disp(['Average Processing Time per Image(s): ', num2str(avgTime(i))]);
end

%% 3) Edge Detection: Sobel
Htime = zeros(1,7); 
avgTime = zeros(1,7);

% Set up Gradient Filters:
Gblur = [1 2 1; 2 4 2; 1 2 1];
Gscale = sum(sum(Gblur));

Hx = [-1 0 1; -2 0 2; -1 0 1]; % X-direction
Hy = [-1 -2 -1; 0 0 0; 1 2 1]; % Y-direction
Scale = 1; 

for i = 1:numClass
    tic 
    numImg = length(CellType{i});
    
    ImgSobel = {};
    HSobel = {};
    
    for j = 1:numImg
      Img = double(Class(i).image{j});
      [Img, Hg] = LinearFilter(Img, Gblur, Gscale);
      [ImgX, HX] = LinearFilter(Img, Hx, Scale);
      [ImgY, HY] = LinearFilter(Img, Hy, Scale);
      ImgSobel{j} = sqrt(abs(ImgX).^2 + abs(ImgY).^2);
      ImgSobel{j} = uint8(ImgSobel{j});
    end
    Class(i).ImgSobel = ImgSobel;

    Htime(i) = toc;
    avgTime(i) = Htime(i) / numImg;
    disp(['Total Processing Time (s): ', num2str(Htime(i))]);
    disp(['Average Processing Time per Image(s): ', num2str(avgTime(i))]);
end

%% 3) Erosion and Dilation Operators
% Convert Imgsobel to Binary
for i = 1:numClass
    tic 
    numImg = length(CellType{i});
    ImgBin = {};
    
    for j = 1:numImg
      ImBin = zeros(W,L);
      Img = Class(i).ImgSobel{j};
      for u = 1:W
          for v = 1:L % find all pixels that are greater than set threshold
                if Img(u,v) >= 50
                    ImBin(u,v) = 1;
                end
            end
      end
    
      ImgBin{j} = ImBin;
          
    end
    Class(i).ImgBin = ImgBin;

    Htime(i) = toc;
    avgTime(i) = Htime(i) / numImg;
    disp(['Total Processing Time (s): ', num2str(Htime(i))]);
    disp(['Average Processing Time per Image(s): ', num2str(avgTime(i))]);

end
%% Dilate Image
B = [0 0; 1 1; -1 1; 1 -1; -1 -1]; % struct element

for i = 1:numClass
    tic 
    numImg = length(CellType{i});
    ImgDil = {};
    
    for j = 1:numImg
      % Dilate Image
      ImgBin = Class(i).ImgBin{j};
      [ix, iy, v] = find(ImgBin == 1);
      newDil = zeros(W,L);
      for b = 1:size(B,1)
          newMask = zeros(W,L);
          for c = 1: length(ix)
            newMask(ix(c) + B(b,1), iy(c) + B(b,2)) = ImgBin(ix(c), iy(c));
          end
          newDil = newDil + newMask; % overlay all tranposed images
      end
      
      ImgDil{j} = newDil; 
      
    end
    
    Class(i).ImgDilate = ImgDil;
    
    Htime(i) = toc;
    avgTime(i) = Htime(i) / numImg;
    disp(['Total Processing Time (s): ', num2str(Htime(i))]);
    disp(['Average Processing Time per Image(s): ', num2str(avgTime(i))]);
end

%% Erode Image
B = [0 1 0; 1 0 1; 0 1 0]; % struct element

% set index based on struct element used
c=floor(size(B,1)/2);
d=floor(size(B,2)/2);

for i = 1:7
    tic 
    numImg = length(CellType{i});
    ImgErode = {};
    for j = 1:numImg
      % Erode Image
      ImgBin = Class(i).ImgBin{j};
      newErode = zeros(W,L);
      for u = 1:W-(2*c)
          for v = 1:L-(2*d)     
              newMask = ImgBin(u:u+(2*c), v:v+(2*d));      
              Ediff = min(min(newMask - B)); 
              if Ediff == 0 
                  % keep if there is no difference between image and B
                 newErode(u,v) = 1;
              end
          end
      end  
      
      ImgErode{j} = newErode;
    end
    
    Class(i).ImgErode = ImgErode; 
    Htime(i) = toc;
    avgTime(i) = Htime(i) / numImg;
    disp(['Total Processing Time (s): ', num2str(Htime(i))]);
    disp(['Average Processing Time per Image(s): ', num2str(avgTime(i))]);
end

%% Image Segmentation: Histogram Thresholding
for i = 1:7    
    tic 
    numImg = length(CellType{i});
    ImgThreshold = {};
    for j = 1:numImg
        Hist = Class(i).hist{j};
        Image = Class(i).image{j};
        
        P = (Hist ./ (W*L))'; % probability of histogram intensities
        Int = 0:255; % vector of pixel intensities
        cdf = cumsum(P); % cumulative density function
        Thresh = 1; % initialize threshold value
        maxV = 0; % initialize max variance value

        for T = 2:255 % run through each intensity value
            obj = Int(1:T); %cells
            back = Int(T+1:end); %background
            Ub(T) = sum(back.*P(T+1:end)) / (1-cdf(T)); %mean cells
            Uo(T) = sum(obj.*P(1:T)) / cdf(T); % mean background
            VarB(T) = cdf(T)*(1-cdf(T))*(Uo(T) - Ub(T)).^2;
                % variance difference between cells and background

            if VarB(T) >= maxV % find the minimum variance difference
                maxV = VarB(T);
                Thresh = T;
            end
        end
        
        % segment image based on the threshold found
        for u = 1:W
            for v = 1:L
                if Image(u,v) <= Thresh % Cell
                    SegImg(u,v) = 1;
                elseif Image(u,v) > Thresh % Background
                    SegImg(u,v) = 0;
                end        
            end
        end
        
        ImgThreshold{j} = SegImg;
    end
    
    Class(i).ImgThreshold = ImgThreshold;
    
    Htime(i) = toc;
    avgTime(i) = Htime(i) / numImg;
    disp(['Total Processing Time (s): ', num2str(Htime(i))]);
    disp(['Average Processing Time per Image(s): ', num2str(avgTime(i))]);
end

%% Image Segmentation: K mean Clustering
k = 5; % number of clusters used
centers = 1:k; % calculate center of clusters

% Set thresholds for the k parameter used (initial clusters and output
% intensities to separate clusters)
thresh = round((W*L) / k); 
indx = 1:thresh:(W*L);
indx(k+1) = W*L;
Ints = linspace(0, 255, k); % new output intensities
Ints = round(Ints);

for i = 1:7
    tic 
    numImg = length(CellType{i});
    ImgKclust = {};

    for j = 1:numImg
        flag = 0; 
        t = 1;
        
        Data = Class(i).image{j};
        Data = double(Data(:));
        for c = 1:k
            Kmean(c) = mean(Data(indx(c:c+1))); % mean given subset of data
        end
        KmeanNew = Kmean;
        while flag == 0 
            % loop through k-means cluster algorithm until cluster means don't change
            Kmean = KmeanNew;
            for u = 1:(W*L)
                for v = 1:k
                    Diff(v) = sqrt((Kmean(v) - Data(u)).^2); 
                        % find euclidean distance between mean and data
                end
                [val, idx(u)] = min(Diff); 
                    % assign to a cluster based on minimum difference      
            end

            for v = 1:k % calculate new cluster means
                cluster = find(idx == v);
                KmeanNew(v) = mean(Data(cluster));
                Comp{v} = cluster;
            end

            if sum(Kmean - KmeanNew) == 0 || t > 10 
                % if no change in mean values, stop looping
                flag = 1;
            end

            t = t + 1; % imposing restriction to prevent infinite looping
        end
        
        % Making the new image based on the clusters found
        NewImage = zeros(W,L);

        for v = 1:k
            idx = Comp{v};
            NewImage(idx) = Ints(v);
        end
        
        ImgKclust{j} = uint8(NewImage);       
    end
    
    Class(i).ImgKcluster = ImgKclust;
    
    Htime(i) = toc;
    avgTime(i) = Htime(i) / numImg;
    disp(['Total Processing Time (s): ', num2str(Htime(i))]);
    disp(['Average Processing Time per Image(s): ', num2str(avgTime(i))]);
end









