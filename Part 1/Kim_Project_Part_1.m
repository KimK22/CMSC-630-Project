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

% Quantization Level
Level = 7;  % Threshold for the # bins wanted

%  Noise Addition Variable Setups
Dnoise    = 0.01;    % noise addition density for salt and pepper filter

Gmean     = 10;      % Gaussian Mean
Gvar      = 50;      % Gaussian Variance

%  Setting up the Spatial Neighborhood
Fweight   = [1 1 1; 1 1 1; 1 1 1]; % Integrated weight vector
Fscale    = sum(Fweight(:));  % Scaling factor for the weights (sum = 1)

Mweight   = [ 1 2 1; 2 3 2; 1 2 1]; % weight vector for median filter


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
    Hsum = zeros(256,1);
    numImg = length(CellType{i});
    for j = 1:numImg
       [ImgEq{j}, Heq{j}] = HistEqual(Class(i).image{j}, Class(i).hist{j});
       Hsum = Hsum + Heq{j}; 
    end
    Class(i).equalImg = ImgEq;
    Class(i).equalH = Heq;
    
    Havg(:,i) = round(Hsum ./ numImg);
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
    ImgQu = {};
    Hq = {};
    Hsum = zeros(1,256);
    MSQE = zeros(1, numImg);
    
    for j = 1:numImg
       [ImgQu{j}, Hq{j}, MSQE(j)] = ImageQuant(Class(i).image{j}, Level);
        Hsum = Hsum + Hq{j}; 
    end
    Class(i).quantImg = ImgQu;
    Class(i).quantH   = Hq;
    Class(i).msqe     = MSQE;
    Havg(:,i) = round(Hsum ./ numImg);
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
    Hsum = zeros(1,256);
    
    for j = 1:numImg
       [ImgSP{j}, Hsp{j}] = SaltPepperNoise(Class(i).image{j}, Dnoise);
        Hsum = Hsum + Hsp{j}; 
    end
    Class(i).ImgSP01 = ImgSP;
    Class(i).Hsp01   = Hsp;
    Havg(:,i) = round(Hsum ./ numImg);
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
    ImgG = {};
    Hgauss = {};
    Hsum = zeros(1,256);
    
    for j = 1:numImg
       [ImgG{j}, Hgauss{j}] = GaussNoise(Class(i).image{j}, Gmean, Gvar);
       Hsum = Hsum + Hgauss{j};
    end
    Class(i).ImgGauss1050 = ImgG;
    Class(i).Hgauss1050   = Hgauss;
    
    Havg(:,i) = round(Hsum ./ numImg);
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
    Hsum = zeros(1,256);
    
    for j = 1:numImg
   [ImgLin{j}, HLin{j}] = LinearFilter(Class(i).image{j}, Fweight, Fscale);
       Hsum = Hsum + HLin{j};
    end
    Class(i).ImgLineargauss = ImgLin;
    Class(i).HLinearguass   = HLin;
    
    Htime(i) = toc;
    Havg(:,i) = round(Hsum ./ numImg);
    avgTime(i) = Htime(i) / numImg;
    disp(['Total Processing Time (s): ', num2str(Htime(i))]);
    disp(['Average Processing Time per Image(s): ', num2str(avgTime(i))]);
end

%% 8) Filtering Operations: Median Filter
Htime = zeros(1,7); 
avgTime = zeros(1,7);

for i = 1:numClass
    tic

    numImg = length(CellType{i});
    ImgMed = {};
    HMed = {};
    Hsum = zeros(1,256);
    
    for j = 1:numImg
       [ImgMed{j}, HMed{j}] = MedFilter(Class(i).image{j}, Mweight);
       Hsum = Hsum + HMed{j};
    end
    Class(i).ImgMedian = ImgMed;
    Class(i).HMedian   = HMed;
    
    Htime(i) = toc;
    Havg(:,i) = round(Hsum ./ numImg);
    avgTime(i) = Htime(i) / numImg;
    disp(['Total Processing Time (s): ', num2str(Htime(i))]);
    disp(['Average Processing Time per Image(s): ', num2str(avgTime(i))]);
end