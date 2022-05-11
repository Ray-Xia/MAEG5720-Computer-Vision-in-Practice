clear, clc, close all

% Preparation
run('vl_setup'); % one-time setup

%% Stitch 2 images
im1 = imread('cu1.JPG');
im2 = imread('cu2.JPG');

% im1 = imread('sceneL.jpg');
% im2 = imread('sceneR.jpg');

im1 = imresize(im1, 0.5); im1 = rgb2gray(im1); im1 = single(im1);
im2 = imresize(im2, 0.5); im2 = rgb2gray(im2); im2 = single(im2);

% Find SIFT Features and Descriptors
[f1,d1] = vl_sift(im1); % k is a 4 (x; y; sigma; theta) x num_features
[f2,d2] = vl_sift(im2); % d is a 128 x num_features descriptor

tic

[matches,scores] = match_descriptor(d1,d2);

[bestH,num_of_inliers] = RANSAC(f1,f2,matches);

stitchedImage= stitch(im1,im2,bestH);

toc

figure('Name','Stiched Image of 2 images');
imshow(stitchedImage);

%% Stitch 3 images
imgLeft = imread('treeL.jpg');
imgCenter = imread('treeC.jpg');
imgRight = imread('treeR.jpg');
tic 
stitch_3_images = main_multiple_images(imgLeft,imgCenter,imgRight); % imgCenter is the reference image
toc 
figure('Name','Stich 3 images');
imshow(stitch_3_images);