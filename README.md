# CMSC-630-Project
Repository that holds all files for CMSC 630 project. 

Part 1: Image processing and filtering
Due:    March 13, 2022

Goal: Design a program that works in a batch setting. Must have a main file that initializes all images given the location of the images, and processes them given the required functions. 

Requirements:
  1) General framework for processing all images in a batch setting with supplying setup initializing file
  2) Noise Addition functions to corrupt images: Salt and Pepper ; Gaussian Noise
  3) Convert Color images to selected single color spectrum
  4) Histogram calculation for each image
  5) averaged histograms of pixel values for each CLASS of images
  6) Selected image quantization technique for user-specified levels
  7) Filter Options: Linear (mask size and weight) ; Median (mask size and weight)
  8) Display Performance Measurements: Processing time of batch and average per image;  MSQE for image quantization levels

Part 2: Image Segmentation
Due:    April 10, 2022

Goal: Design a program that can implement different image segmentation techniques. Must work in a batch setting given the location of the images. 

Requirements: 
  1) Implement one selected edge detection algorithm
  2) Implement dilation adn erosion operators
  3) Impement two segmentation techniques
        - Histogram Thresholding: single threshold that divides the image into 2 segments: foreground (cells) and background (everything else)
        - K-means Clustering: examine effect of different values of k parameter on the segmentation
