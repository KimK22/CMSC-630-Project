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

ImgOrig = imread(mainFolder(3).name);
ImgRed = ImgOrig(:,:,1);
ImgGreen = ImgOrig(:,:,2);
ImgBlue = ImgOrig(:,:,3);
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

%% Setting up Variables
%  Choosing Color Channel
Color = 1; % Red = 1; Green = 2; Blue = 3

%  Setting up the Spatial Neighborhood
Fsize = 3; % size of the matrix : assume nxn matrix size
Fweight = [ 1 1 1; 5 5 5; 1 1 1]; % weight of spatial neighborhood
Fscale = sum(Fweight(:)); % Scaling factor for the weights (sum = 1)

%% Histogram Calculations
% Histogram calculation for each original grayscale image
Img = ImgRed;
h = zeros(256,1);
maxInt = 255;
[W, L] = size(Img);
PixTotal = W*L; % total number of pixels

for i = 1:maxInt+1
    h(i) = sum(sum(ImgQ == i));
    pdf(i) = h(i) / (W*L); % probability density (# pixels of intensity H / total # pixels)
end
subplot(2,2,1); imshow(ImgQ); subplot(2,2,2);
bar(0:255,h); % plot histogram

%% Histogram Equalization
% Goal: Compress based on peak values in the histogram 
ImgEq = ImgRed; % copy of original image
total = 0;
for i = 1:maxInt+1
    total = total + h(i);
    cumtotal(i) = total; 
    cdf(i) = cumtotal(i) / PixTotal;
    Out(i) = round(cdf(i)*maxInt);
end

for i = 1:W
    for j = 1:L
        ImgEq(i,j) = Out(Img(i,j) + 1);
    end
end

subplot(2,2,3); imshow(ImgEq); subplot(2,2,4);
bar(0:255, Out);

%% Histogram Quantization
% Goal: Reduce number of bins in histogram / reduce values a pixel can be
% in the image.
Levels = 20;
maxInt = 255;
spaces = floor(maxInt / Levels);
Thresholds = 0:spaces:maxInt;
for val = 2:Levels
    vals = [Thresholds(val-1):Thresholds(val)];
    summ = sum(vals);
    total = length(vals);
    PixInt(val) = summ/total;
end
Img = ImgRed;
ImgQ = Img;
[W, L] = size(Img);
for i = 1:W
    for j = 1:L
        for val = 1:Levels-1
            if Img(i,j) >= Thresholds(val) && Img(i,j) < Thresholds(val+1)
                ImgQ(i,j) = PixInt(val);
            end
        end
    end
end
subplot(2,1,1); imshow(Img);
subplot(2,1,2); imshow(ImgQ);


%% Filtering Options
ImgLin = Img; 

% 1) Linear Filter
for u = 2:W-1
    for v = 2:L-1
        Psum = 0; % sum of average
        for i = -1:1
            for j = -1:1
                Psum = Psum + Img(u+i, v+j)*Fweight(i+2, j+2)/Fscale;
            end
        end
        ImgLin(u,v) = round(Psum);
    end
end

imshow(ImgLin);

%% 2) Median Filter
ImgMed = Img; % copy original Image
MedFilt = [1 4 2; 3 5 3; 0 0 0];

for u = 2:W-1
    for v = 2:L-1
        inputs = [];
        for i = -1:1
            for j = -1:1
               for val = 1:MedFilt(i+2, j+2)   
                   inputs = [inputs, Img(u+i,v+j)];
               end
            end
        end
        ImgMed(u,v) = median(inputs);       
    end
end

imshow(ImgMed);

%% Noise Additions
% 1) Salt and Pepper noise of user-specified strength
ImgSalt = Img;

Dnoise = 0.005; % percent of pixels affected by noise
Noise = rand(W, L);
for i = 1:W
    for j = 1:L
        if Noise(i,j) > 0 && Noise(i,j) < Dnoise / 2
            ImgSalt(i,j) = 0;
        elseif Noise(i,j) > Dnoise/2 && Noise(i,j) < Dnoise
            ImgSalt(i,j) = maxInt;
        end         
    end
end

imshow(ImgSalt);

%% 2) Gaussian Noise of user-specified parameters
ImgGauss = Img;
Gmean = 1;  % Gaussian Mean
Gvar = 0.1; % Gaussian Variance 

for i = 1:W
    for j = 1:L
        Gnoise = Gmean + randn(1)*sqrt(Gvar);
        ImgGauss(i,j) = Img(i,j) + Gnoise; 
    end
end

imshow(ImgGauss);

%% Processing Times
tic
time = toc; % stores end processing time (sec)
% Can turn into vector to average later on

%% MSQE for image quantization levels











