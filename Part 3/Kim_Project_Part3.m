% Name: Kristin Kim
% Course: CMSC 630
% Outline: Project Part 3 Feature Extraction + Image Classification

% Goal: 
%   Use cell images as training data to classify cell type. 

% Requirements: 
%     1) From segmented cell images (Choose any technique) extract at least
%     4 distinct features and assign class label based on cell type from
%     documentation (7 in total).
%     2) Save new dataset as a matrix with columns = feature in .csv format
%     3) Implement a k-NN classifier with Euclidean distance
%           A) Perform classification and report accuracy 
%           B) Evaluate the performance of parameter k on classification 
%               accuracy (run independent experiments with at least 5 
%               different k values and compare results).
%     4) Implement 10-fold cross-validation
%           A) Perform classification of cells. Report classification 
%               accuracy (averaged among all 10 fold of CV).
%     5) Present details of features you extracted nad implementation
%     specific knn adn obtained results in a report.

% Note: Super045.BMP = corrupted file, so I removed it

%% Cleanup Section
clc; 
close all; 
clear all;

%% 1) Batch Image Initialization
Folder = 'Cancerous cell smears'; % Location of images
addpath(Folder);
mainFolder = dir(Folder); % Stores all image names in struct format
mainFolder(1:3) = [];
numImg = length(mainFolder);
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

% Setup Variables
%  Integrate categories into same cell
CellType  = {cyl, inter, let, mod, para, super, svar};
CellClass = [1 2 3 4 5 6 7]; % classification of cell
CellName  = {'cyl', 'inter', 'let', 'mod', 'para', 'super', 'svar'};
%  Choosing Color Channel
Color     = 1;      % Red = 1; Green = 2; Blue = 3
% Size of Image
[W,L] = size(cyl{1}(:,:,1)); % general size of images (Assumed constant)
% Set up Gradient Filters:
Gblur = [1 2 1; 2 4 2; 1 2 1];
Gscale = sum(sum(Gblur));
% Binary Image Erosion and Dilation struct elements
E = [0 1 0; 1 1 1; 0 1 0]; % erode image
D = [0 0; 1 1; -1 1; 1 -1; -1 -1]; % struct element

Class = struct();
numClass = length(CellType);

% Set up Cell Image structure
for i = 1:numClass
    Class(i).name  = CellName{i};
    CellColor = {};
    H = {};
    numImg = length(CellType{i});
    for j = 1:numImg
        Ind = cell2mat(CellType{i}(j));
        CellColor{j} = Ind(:,:,Color);
        H{j} = Histogram(CellColor{j});
    end
    
    Class(i).image = CellColor;
    Class(i).hist  = H;
end

%% 2) Image Segmentation and Feature Extraction
CellFeat = zeros(499, 6);
count = 1;

for Nclass = 4% :numClass
    numImg = length(CellType{Nclass});
    CellImg = {};
    for nImg = 1:numImg
        Image = Class(Nclass).image{nImg};
        [GImg] = KmeanCluster(Image, 5);

        Imbin = zeros(W,L);
        for i = 1:W
            for j = 1:L
                if GImg(i,j) < 170 % Identify cells in the image
                    Imbin(i,j) = 1;
                end
            end
        end
        
        Eimg = ErodeImg(Imbin, E);
        Eimg = ErodeImg(Eimg, E);
        
        Dimg = DilateImg(Eimg, D);
        Dimg = DilateImg(Dimg, D);

        Edge = Dimg - Eimg; % Cell perimeters
        Edge = DilateImg(Edge, D);
        Edge = ErodeImg(Edge, E);
        
        TotalArea = sum(sum(Imbin)); % find total cell area in the image
        
        if (TotalArea / (W*L)) > 0.4 % If cells can't be isolated
            Cbody = Imbin;
            Coutline = Edge;
        else
            [Cbody, Coutline] = CellIsolation(Edge); %single cell isolation
        end
        
        % FeatExtraction; CellArea, Perimeter, AvgInt, MaxInt, NucleiArea]
        [ImgFeat, maskImg, maskImgC] = CellFeatExtract(Cbody, Coutline, Image, GImg);
        
        % Save feats in matrix
        CellFeat(count, :) = [ImgFeat Nclass];
        count = count + 1;
        
        CellImg{nImg} = maskImg;
        CellImgC{nImg} = maskImgC;
        TotArea{nImg} = Cbody;
        TotPerim{nImg} = Coutline;
    end
        Class(Nclass).SegmentedCell = CellImg;
        Class(Nclass).SegmentedCellC = CellImgC;
        Class(Nclass).CellArea = TotArea;
        Class(Nclass).CellPerim = TotPerim;
end

%% Save features in excel file
filename = 'CellFeatureData.xlsx';
writematrix(CellFeat, filename);

%% 3) Image Classification based on extracted features

% Import data as a matrix with features adn class labels from a .csv file
CellFeat = xlsread('CellFeatureData.xlsx');
% Values: [Area, Perimeter, Avg Intensity, Max Intensity, Nucleus Area,
% True Class Value]

% Create fold validation groups
Group = {};
TrueClass = {};
ind = [1, 51, 101, 201, 301, 351, 401];
count = [5, 5, 10, 10, 5, 5, 10];
for i = 1:10
    Group{i} = [];
    TrueClass{i} = [];
    for j = 1:length(ind)
       Group{i} = [Group{i}; CellFeat(ind(j):ind(j)+(count(j)-1), 1:5)];
       TrueClass{i}=[TrueClass{i}; CellFeat(ind(j):ind(j)+(count(j)-1),6)];
       ind(j) = ind(j) + count(j);
    end
end

% 10-Fold Cross Validation: 1 data group = testing data, the other 9 will
% be the training data that you're comparing to. Conduct K-nn for all
% combinations and observe results.

% Set up the 10-fold cross validation test and training groups
TestData = {};
TrainData = {};
TrainClass = {};
for fold = 1:10
    TestData{fold} = [];
    TrainData{fold} = [];
    TrainClass{fold} = [];
    for i = 1:10
        if i == fold
            TestData{fold} = Group{i};           
        else
            TrainData{fold} = [TrainData{fold}; Group{i}];
            TrainClass{fold} = [TrainClass{fold}; TrueClass{i}];
        end
    end
end

% K-Nearest Neighbor Euclidean Distance Detection
k = [2, 4, 5, 10, 15, 20]; 
acc = [];
avg_acc = zeros(1, length(k));

for Kval = 1:length(k) % loop through different k-parameter values
    for i = 1:10
        acc(i) = 0; % accuracy per fold
        TestD = TestData{i};
        TrainD = TrainData{i};
        TrainC = TrainClass{i};
        TrueC = TrueClass{i};
        PredictC = zeros(1,50);
        for n = 1:length(TestD)
            dist= sqrt(sum((TrainD - TestD(n,:)).^2')); 
                % distance between Test and Train
            [dist2, ind] = sort(dist, 'ascend'); % find minimum distances
            TrainMaj = TrainC(ind); % resorting Training class values
            class_pos = TrainMaj(1:k(Kval)); % find k minimum values
            PredictC(n) = mode(class_pos); % estimate mode(Test class)

            if (PredictC(n) == TrueC(n))
                acc(i) = acc(i) + 1;
            end
        end
        acc(i) = acc(i)/length(TestD); % accuracy of the fold validation
    end

    avg_acc(Kval) = mean(acc); % total accuracy for k-value
end