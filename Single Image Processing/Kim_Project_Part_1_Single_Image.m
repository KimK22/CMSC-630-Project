% Name: Kristin Kim
% Course: CMSC 630
% Outline: Project Part 1

% Goals: Call 500 real-life cancerous smear images and edit them using
% different functionalities. This will be the main setup file that has the
% folder/ images imported and runs the functions defined. 

% Will have separate .m files that code for each function that will be run
% through the main code file. 
%% Note: Super045.BMP = corrupted file, so I removed it

%% Cleanup Section
clc;
close all;
clear all;

%% 1) General framework for processing all images in a batch setting with supplying set-up initializing file
% Open image folder and add to search path
Folder = 'Cancerous cell smears';
addpath(Folder);
mainFolder = dir(Folder); % Stores all image names in struct format

%% Separate different classes into variables
% for i = 3:length(mainFolder)
%     if mainFolder(i).name(1:3) == 'cyl'
%         cyl{i-2} = imread(mainFolder(i).name);
%     elseif mainFolder(i).name(1:3) == 'int'
%         inter{i-52} = imread(mainFolder(i).name);
%     elseif mainFolder(i).name(1:3) == 'let'
%         let{i - 102} = imread(mainFolder(i).name);
%     elseif mainFolder(i).name(1:3) == 'mod'
%         mod{i - 202} = imread(mainFolder(i).name);
%     elseif mainFolder(i).name(1:3) == 'par'
%         para{i - 302} = imread(mainFolder(i).name);
%     elseif mainFolder(i).name(1:3) == 'sup'
%         super{i - 352} = imread(mainFolder(i).name);
%     elseif mainFolder(i).name(1:3) == 'sva'
%         svar{i - 401} = imread(mainFolder(i).name);
%     end
% end

%% Setup Variables
%  Integrate categories into same cell
% Classes   = {cyl, inter, let, mod, para, super, svar};

%  Choosing Color Channel
Color     = 2;      % Red = 1; Green = 2; Blue = 3

%  Setting up the Spatial Neighborhood
Fsize     = 3;      % size of the matrix : assume nxn matrix size
Fweight   = [ 0 1 0; 1 4 1; 0 1 0]; % weight of spatial neighborhood
Fscale    = sum(Fweight(:)); % Scaling factor for the weights (sum = 1)
Ftype     = 1;      % Filter type: 1) Linear, 2) Median

Mweight   = [ 0 4 0; 3 0 3; 0 4 0];
Mscale    = 1;

%  Histogram Variable Setups
SumHist   = 0;      % Sum of histogram values for each class. 
Level     = 7;      % Quantization Level threshold
%  Noise Addition Variable Setups
Dnoise    = 0.01;   % noise addition for salt and pepper filters
Gmean     = 10;      % Gaussian Mean
Gvar      = 50;    % Gaussian Variance

%% Running the Main Code

% Original Image: Color Channel Chosen
Image = mainFolder(200).name;
Image = imread(Image);
subplot(5,2,1); imshow(Image(:,:,Color));
[H, pdf] = Histogram(Image, Color); subplot(5,2,2); bar(0:255, H);

% Linear Filtering
[LinImg, Hlin] = ImgFilter(Image, Color, Fweight, Fscale, Ftype);
subplot(5,2,3); imshow(LinImg);
subplot(5,2,4); bar(0:255, Hlin);

% Median Filtering
[MedImg, Hmed] = ImgFilter(Image, Color, Mweight, Mscale, 2);
subplot(5,2,5); imshow(MedImg);
subplot(5,2,6); bar(0:255, Hmed);

% Salt and Pepper Noise Addition
[SPImg, Hsp] = SaltPepperNoise(Image, Color, Dnoise);
subplot(5,2,7); imshow(SPImg);
subplot(5,2,8); bar(0:255, Hsp);

% Gaussian Noise Addition
[GNImage, Hg] = GaussNoise(Image, Color, Gmean, Gvar);
subplot(5,2,9); imshow(GNImage);
subplot(5,2,10); bar(0:255, Hg);

%% Histogram Calculations
figure(); 
% Original Image
subplot(4,2,1); imshow(Image(:,:,Color));
[H, pdf] = Histogram(Image, Color); subplot(4,2,2); bar(0:255, H);

% Image Quantization
[ImgQ, Hq] = ImageQuant(Image, Color, Level);
subplot(4,2,3); imshow(ImgQ);
subplot(4,2,4); bar(0:255, Hq);

% Histogram Equalization
[ImgEq, Heq] = HistEqual(Image, Color, H);
subplot(4,2,7); imshow(ImgEq);
subplot(4,2,8); bar(0:255, Heq);

%% Trying to use each against each other
figure();
subplot(5,2,1); imshow(Image(:,:,Color));
[H, pdf] = Histogram(Image, Color); subplot(5,2,2); bar(0:255, H);

% 1) ADd noise to the image
Dnoise = 0.01;
[NoisyImage, Hnoise] = SaltPepperNoise(Image, Color, Dnoise);
%NoisyImage = GaussNoise(Image, Color, Gmean, Gvar); 
subplot(5,2,3); imshow(NoisyImage);
subplot(5,2,4); bar(0:255, Hnoise);

% 2) Remove Noise using Linear Filter
[FSPImg, Hf] = ImgFilter(NoisyImage, 1, Fweight, Fscale, 1);
subplot(5,2,5); imshow(FSPImg);
subplot(5,2,6); bar(0:255, Hf);

% 3) Using a Median Filter
[MSPImg, Hm] = ImgFilter(NoisyImage, 1, Mweight, Mscale, 2); 
subplot(5,2,7); imshow(MSPImg);
subplot(5,2,8); bar(0:255, Hm);


