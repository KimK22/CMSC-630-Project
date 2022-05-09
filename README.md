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

Part 3: Feature extraction and image classification
Due:   May 12, 2022

Goal: design a program that isolates individual cells in each image. Perform feature extraction that can be used to classify cell types.

Requirements:
  1) From segmented cell images, extract at least 4 distinctive features and assign class label according to cell type from documentation. 
  2) Save new dataset as a matrix with columns representing features, and last as class labels. Save as .csv format
  3) Implement a k-NN classifier with Euclidean distance
  4) Implement 10-fold cross-validation
  5) Perform classification of cells using 10-fold cross-validation adn k-NN classifier, Report classification accuracy (Averaged among 10 folds of cross validation)
  6) Evaluate performance of parameters using 5 different values of K
