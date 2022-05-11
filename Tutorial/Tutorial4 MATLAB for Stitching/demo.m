clear, clc, close all
% Preparation
run('vl_setup'); % run once

im1 = imread('cu1.JPG');
im2 = imread('cu2.JPG');

im1 = imresize(im1, 0.5); im1 = rgb2gray(im1); im1 = single(im1);
im2 = imresize(im2, 0.5); im2 = rgb2gray(im2); im2 = single(im2);

% Find SIFT Features and Descriptors
[k1,d1] = vl_sift(im1); % k is a 4 (x,y,sigma, theta) x num_features
[k2,d2] = vl_sift(im2); % d is a 128 x num_features

%% Hints for feature match: calculate the SSD between features
% match the first feature in img1
d1_selected = d1(:, 1);
temp = zeros(2,size(k2,2)); thresh = 1.5; k = 0; scores = zeros(1,size(k1,2)); matches = zeros(2,size(k1,2));
for j = 1 : size(k2,2)
    temp(1 , j) = j;
    temp(2 , j) = norm(double(d1(:,1)-d2(:,j)),2);
end
% sort on the 3rd row
sorted = transpose(sortrows(transpose(temp),2));
min1 = sorted(:,1);
min2 = sorted(:,2);
min2(2,:) / min1(2,:) >= thresh

%% Hints for robust H estimation with RANSAC
% Take the initial match as example here
load('matches.mat') % need to be done by your self.
numMatches = size(matches, 2); numPts = 4;
randSample = randperm(numMatches, numPts);
d1Ind_pt1 = matches(1, randSample(1));
d2Ind_pt1 = matches(2, randSample(1));
d1Ind_pt2 = matches(1, randSample(2));
d2Ind_pt2 = matches(2, randSample(2));
d1Ind_pt3 = matches(1, randSample(3));
d2Ind_pt3 = matches(2, randSample(3));
d1Ind_pt4 = matches(1, randSample(4));
d2Ind_pt4 = matches(2, randSample(4));
H = FindHomography(k1(1, d1Ind_pt1), k1(2, d1Ind_pt1), k1(1, d1Ind_pt2), k1(2, d1Ind_pt2), k1(1, d1Ind_pt3), k1(2, d1Ind_pt3), k1(1, d1Ind_pt4), k1(2, d1Ind_pt4), ...
    k2(1, d2Ind_pt1), k2(2, d2Ind_pt1), k2(1, d2Ind_pt2), k2(2, d2Ind_pt2), k2(1, d2Ind_pt3), k2(2, d2Ind_pt3), k2(1, d2Ind_pt4), k2(2, d2Ind_pt4));

%% warp pts
matchIndex = 10; % need to loop all matched pts
d1Ind = matches(1, matchIndex);
d2Ind  = matches(2, matchIndex);

X1 = k1(1, d1Ind); Y1 = k1(2, d1Ind);
X2  = k2(1, d2Ind); Y2  = k2(2, d2Ind);

P1 = [X1; Y1; 1];
P2  = [X2; Y2; 1];

P11 = H*P1;
P11 = P11 ./ P11(3);
err = norm((P11 - P2),2); % compare with err <= pixelError to get the inliners / use the largest inliners to get the best H

%% stitch image using the best H
load('bestHomography.mat') % need to be done by your self.
stitchedImage = im1;
stitchedImage = padarray(stitchedImage, [0 size(im2, 2)], 0, 'post');
stitchedImage = padarray(stitchedImage, [size(im2, 1) 0], 0, 'both');

i = 1400; j = 600;
p2 = bestHomography * [i; j-size(im2, 1); 1];
p2 = p2 ./ p2(3);

x2 = round(p2(1));
y2 = round(p2(2));

if x2 > 0 && x2 <= size(im2, 2) && y2 > 0 && y2 <= size(im2, 1)
    if (im2(y2, x2)>stitchedImage(j, i))         
        stitchedImage(j, i) = im2(y2, x2);
    end

end
imshow(stitchedImage,[])

% remove the extra boundary (which is 0)
% use MATLAB function: find 
% then imcrop