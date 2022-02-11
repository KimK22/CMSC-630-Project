% Name: Kristin Kim
% Course: CMSC 630
% Outline: Project Part 1

% Goals: Call 500 real-life cancerous smear images and edit them using
% different functionalities. This will be the main setup file that has the
% folder/ images imported and runs the functions defined. 

% Initial Coding: Use only 1 image input

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

Img0 = imread(mainFolder(3).name);
figure(1); subplot(2,2,1); imshow(Img0);

% 2) Convert color images to grayscale
Img = Img0(:,:,1)*0.3 + Img0(:,:,2)*0.59 + Img0(:,:,3)*0.11;
subplot(2,2,2); imshow(Img);

% 3) Histogram calculation for each original grayscale image
h = zeros(256,1);
for I = 0:255
    h(I + 1) = sum(sum(Img(:,:,1) == I));
end
bar(0:255,h); % plot histogram

%% 4) Noise addition functions to corrupt images with: Salt and pepper , Gaussian Noise

%% 5) Averaged histograms of pixel values for each class of images (7 total)

%% 6) Histogram equalization for each image

%% 7) Selected image quantization technique (1 per class)

%% 8) Filtering operations (Linear, Median)
tic 
% Linear Filter given user input
Fgaussian = [0 1 0; 2 5 2; 0 1 0]; % gaussian weight
scale = sum(Fgaussian(:)); % Scaling factor
% [Img2] = ImgFilter(Img, Fgaussian, scale);

val = size(Img);
w = val(1); % width of Image
h = val(2); % height of Image

Img2 = Img; % copy original Image

for u = 2:w-1
    for v = 2:h-1
        Psum = 0; % sum of average
        for i = -1:1
            for j = -1:1
                Psum = Psum + Img(u+i, v+j)*Fgaussian(i+2, j+2)/scale;
            end
        end
        Img2(u,v) = round(Psum);
    end
end

subplot(2,2,3); imshow(Img2);

toc
%%  Median Filter
tic

Img3 = Img; % copy original Image
for u = 2:w-1
    for v = 2:h-1
        inputs = [];
        for i = -1:1
            for j = -1:1
               inputs = [inputs; Img(u+i, v+j)];               
            end
        end
        Img3(u,v) = median(inputs);       
    end
end

subplot(2,2,4); imshow(Img3);

toc
%% 9) Display Performance Measurements:
% Processing time for entire batch per procedure
tic
pause(1)
toc
% Average processing time per image for each
%       tic toc for whole process in each batch, then divide by the total
%       number of images in the batch
% MSQE for image quantization levels

