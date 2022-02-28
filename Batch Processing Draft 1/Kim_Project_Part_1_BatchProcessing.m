% Name: Kristin Kim
% Course: CMSC 630
% Outline: Project Part 1

% Goal: 
%     Call 500 real-life cancerous smear images and edit them using
%     different functionalities. This will be the main setup file that has 
%     the folder/ images imported and runs the functions defined. 

% Requirements: 
%     1) Batch setting image processing given folder location
%     2) Noise addition functions: Salt and Pepper ; Gaussian Noise
%     3) Convert color images to single color spectrum
%     4) Histogram calculation and average per class
%     5) Histogram equalization for EACH image
%     6) Image Quantization 
%     7) Filtering: Linear, Median
%     8) Display Performance Measurements: 
%             Processing time (per procedure and average per image)
%             MSQE for image quantization levels

% Note: Super045.BMP = corrupted file, so I removed it

%% Cleanup Section
clc; 
close all; 
clear all;

%% 1) Batch Image Initialization
Folder = 'Cancerous cell smears'; % Location of images
addpath(Folder);
mainFolder = dir(Folder); % Stores all image names in struct format

% Separate different classes into variables
for j = 3:length(mainFolder)
    if mainFolder(j).name(1:3) == 'cyl'
        cyl{j-2} = imread(mainFolder(j).name);
    elseif mainFolder(j).name(1:3) == 'int'
        inter{j-52} = imread(mainFolder(j).name);
    elseif mainFolder(j).name(1:3) == 'let'
        let{j - 102} = imread(mainFolder(j).name);
    elseif mainFolder(j).name(1:3) == 'mod'
        mod{j - 202} = imread(mainFolder(j).name);
    elseif mainFolder(j).name(1:3) == 'par'
        para{j - 302} = imread(mainFolder(j).name);
    elseif mainFolder(j).name(1:3) == 'sup'
        super{j - 352} = imread(mainFolder(j).name);
    elseif mainFolder(j).name(1:3) == 'sva'
        svar{j - 401} = imread(mainFolder(j).name);
    end
end

%% Setup Variables
%  Integrate categories into same cell
CellType  = {cyl, inter, let, mod, para, super, svar};
CellName  = {'cyl', 'inter', 'let', 'mod', 'para', 'super', 'svar'};
%  Choosing Color Channel
Color     = 1;      % Red = 1; Green = 2; Blue = 3

% Quantization Level
Level = 7;  % Threshold for the # bins wanted

%  Noise Addition Variable Setups
Dnoise    = 0.01;    % noise addition density for salt and pepper filter
Gmean     = 10;      % Gaussian Mean
Gvar      = 50;      % Gaussian Variance

%  Setting up the Spatial Neighborhood
Fsize     = 3;           % size of the matrix : assume nxn matrix size
Wx        = [ 0 1 0];    % weight of spatial neighborhood in x direction
Wy        = [ 1 ; 3; 1]; % weight of spatial neighborhood in y direction
Wt        = Wx + Wy;     % Combined weight matrix

Fweight   = [1 2 1; 3 4 3; 1 2 1]; % Integrated weight vector


Mweight   = [ 0 4 0; 3 0 3; 0 4 0];
Mscale    = 1;

Wscale    = sum(Fweight(:));  % Scaling factor for the weights (sum = 1)
Wtype     = 1;           % Filter type: 1) Linear, 2) Median

%% 2) Histogram Calculations
% Create Class structure with all images, color conversion
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

%% 3) Histogram Equalization
Htime = zeros(1,7); 
avgTime = zeros(1,7);

for i = 1:numClass
    tic
    ImgEq = {};
    Heq = {};
    numImg = length(CellType{i});
    for j = 1:numImg
       [ImgEq{j}, Heq{j}] = HistEqual(Class(i).image{j}, Class(i).hist{j});
    end
    Class(i).equalImg = ImgEq;
    Class(i).equalH = Heq;
    Htime(i) = toc;
    avgTime(i) = Htime(i) / numImg;
    disp(['Total Processing Time (s): ', num2str(Htime(i))]);
    disp(['Average Processing Time per Image(s): ', num2str(avgTime(i))]);
end

%% 4) Image Quantization
Htime = zeros(1,7); 
avgTime = zeros(1,7);

for i = 1:numClass
    tic

    numImg = length(CellType{i});
    ImgSP = {};
    Hsp = {};
    MSQE = zeros(1, numImg);
    
    for j = 1:numImg
       [ImgSP{j}, Hsp{j}, MSQE(j)] = ImageQuant(Class(i).image{j}, Level);
    end
    Class(i).quantImg = ImgSP;
    Class(i).quantH   = Hsp;
    Class(i).msqe     = MSQE;
    
    Htime(i) = toc;
    avgTime(i) = Htime(i) / numImg;
    disp(['Total Processing Time (s): ', num2str(Htime(i))]);
    disp(['Average Processing Time per Image(s): ', num2str(avgTime(i))]);
end

%% 5) Image Noise Addition: Salt and Pepper Noise
Htime = zeros(1,7); 
avgTime = zeros(1,7);

for i = 1:numClass
    tic

    numImg = length(CellType{i});
    ImgSP = {};
    Hsp = {};
    
    for j = 1:numImg
       [ImgSP{j}, Hsp{j}] = SaltPepperNoise(Class(i).image{j}, Dnoise);
    end
    Class(i).ImgSP = ImgSP;
    Class(i).Hsp   = Hsp;
    
    Htime(i) = toc;
    avgTime(i) = Htime(i) / numImg;
    disp(['Total Processing Time (s): ', num2str(Htime(i))]);
    disp(['Average Processing Time per Image(s): ', num2str(avgTime(i))]);
end
%% 6) Image Noise Addition: Gaussian Noise
Htime = zeros(1,7); 
avgTime = zeros(1,7);

for i = 1:numClass
    tic

    numImg = length(CellType{i});
    ImgSP = {};
    Hsp = {};
    
    for j = 1:numImg
       [ImgG{j}, Hlin{j}] = GaussNoise(Class(i).image{j}, Gmean, Gvar);
    end
    Class(i).ImgGauss = ImgG;
    Class(i).Hgauss   = Hlin;
    
    Htime(i) = toc;
    avgTime(i) = Htime(i) / numImg;
    disp(['Total Processing Time (s): ', num2str(Htime(i))]);
    disp(['Average Processing Time per Image(s): ', num2str(avgTime(i))]);
end
%% 7) Filtering Operations: Linear Filter
Htime = zeros(1,7); 
avgTime = zeros(1,7);

for i = 1:numClass
    tic

    numImg = length(CellType{i});
    ImgLin = {};
    Hlin = {};
    
    for j = 1:numImg
       [ImgLin{j}, HLin{j}] = ImgFilter(Class(i).image{j}, Fweight, Wscale, 1);
    end
    Class(i).ImgLinear = ImgLin;
    Class(i).HLinear   = HLin;
    
    Htime(i) = toc;
    avgTime(i) = Htime(i) / numImg;
    disp(['Total Processing Time (s): ', num2str(Htime(i))]);
    disp(['Average Processing Time per Image(s): ', num2str(avgTime(i))]);
end


%% 8) Filtering Operations: Median Filter

% This runs really slow: How can I make it go faster?
Htime = zeros(1,7); 
avgTime = zeros(1,7);

for i = 1:numClass
    tic

    numImg = length(CellType{i});
    ImgMed = {};
    HMed = {};
    
    for j = 1:numImg
       [ImgMed{j}, HMed{j}] = ImgFilter(Class(i).image{j}, Mweight, Mscale, 2);
    end
    Class(i).ImgMedian = ImgMed;
    Class(i).HMedian   = HMed;
    
    Htime(i) = toc;
    avgTime(i) = Htime(i) / numImg;
    disp(['Total Processing Time (s): ', num2str(Htime(i))]);
    disp(['Average Processing Time per Image(s): ', num2str(avgTime(i))]);
end
















