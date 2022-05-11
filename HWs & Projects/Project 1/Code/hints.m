%%%%%%%%%%%%%%%%%%%%%
clc;clear;
currentFolder = pwd;
path_to_images = [currentFolder, '\CarSequence'];
% load the template
fname = sprintf('%s\\car_template.jpg',path_to_images);
img1_rgb = imread(fname);
img1_gray = rgb2gray(img1_rgb);
window = [0,9,9,0]; %x1 y1 x2 y2, the coordinates of top-left and bottom-right corners
%%
%%%%%%%%%%%%%%%%%%%% use Nomarlized Cross-Correlation(NCC) to deal with first frame
frame = 1;
fname = sprintf('%s\\frame00%d.jpg',path_to_images,frame+302);
img2_rgb = imread(fname);
img2_gray = rgb2gray(img2_rgb);

img1_gray = im2double(img1_gray);
img2_gray = im2double(img2_gray);
[row1, col1] = size(img1_gray);
[row2, col2] = size(img2_gray);
% smooth use gaussian kernel --sigma
sigma = 10;
gaussFilt = fspecial('gaussian',10,sigma);  % SEE 'help fspecial' 
img1_gray = imfilter(img1_gray, gaussFilt, 'symmetric'); 
img2_gray = imfilter(img2_gray, gaussFilt, 'symmetric'); 
% Initial bbox info using NCC
I_NCC=normxcorr2(img1_gray,img2_gray);
[ypeak,xpeak]=find(I_NCC==max(I_NCC(:))); % bottom right pt: y2, x2
window = [xpeak-col1+1,ypeak-row1+1,xpeak,ypeak]; % x1,y1,x2,y2
imshow(img2_rgb);
rectangle('Position',[window(1),window(2),col1,row1],'edgecolor','y');
%%
%%%%%%%%%%%%%%%%% or you can use feature-match to deal with the first frame
% frame = 1;
% fname = sprintf('%s\\frame00%d.jpg',path_to_images,frame+302);
% img2_rgb = imread(fname);
% img2_gray = rgb2gray(img2_rgb);
% img2_gray = im2double(img2_gray);
% img2_gray = imfilter(img2_gray, gaussFilt, 'symmetric'); 
T_cps = detectHarrisFeatures(img1_gray, 'MinQuality', 0);
I_cps = detectHarrisFeatures(img2_gray, 'MinQuality', 0);
[features1, valid_points1] = extractFeatures(img1_gray, T_cps);
[features2, valid_points2] = extractFeatures(img2_gray, I_cps);

% 2. match Features
indexPairs = matchFeatures(features1, features2);
matched_points1 = valid_points1(indexPairs(:, 1), :);
matched_points2 = valid_points2(indexPairs(:, 2), :);

% 3. get initial window
p = mean(matched_points2.Location - matched_points1.Location); % top left: x1, y1
p = ceil(p);
window = [p(1), p(2), p(1)+col1-1, p(2)+row1-1];% x1,y1,x2,y2
imshow(img2_rgb);
rectangle('Position',[window(1),window(2),col1,row1],'edgecolor','y');
%%
%%%%%%%%%%%%%%%%%%%% use Lucas Kanade to deal with following frames
% for example we choose the second frame here
window = [330,   213,   417,   261];
frame = 2;
fname = sprintf('%s\\frame00%d.jpg',path_to_images,frame+302);
img2_rgb = imread(fname);
img2_gray = rgb2gray(img2_rgb);
% smooth use gaussian kernel --sigma
gaussFilt = fspecial('gaussian',10,sigma);  % SEE 'help fspecial' 
img2_gray = imfilter(img2_gray, gaussFilt, 'symmetric'); 
img2_gray = im2double(img2_gray);

p = [0; 0];
step = 0;
threshold = 0.01;
max_iteration = 100;

[Tx, Ty] = meshgrid((1:col1)+window(1)-1, (1:row1)+window(2)-1);
for i=1:max_iteration
    % 1. Warp I with W
    warped = interp2(img2_gray, Tx + p(1), Ty + p(2), 'Bilinear');
    warped(isnan(warped)) = 0;
    % 2. Subtract I from T
    I_error = img1_gray(:) - warped(:);
    % 3. Compute Gradient
%     [Gx, Gy] = imgradientxy(warped, 'central');
    [Gx, Gy] = imgradientxy(img1_gray, 'central');
    Ix = Gx(:);
    Iy = Gy(:);
    % 4. Evaluate the jacobian
    W = [1, 0; 0, 1];
    % 5. Compute steepest decent
    A = [Ix, Iy] * W;
    % 6. Compute inverse Hessian
    invH = pinv(A' * A);
    % 7. Multiply steepest descend with error and 8. compute dp
    dp =  invH * A' * I_error;
    % 9. Update the parameters p <- p + delta_p
    p = p + dp;
    if norm(dp) < threshold
        break;
    end
    step = step + 1;
end
window = window + [p(1), p(2), p(1), p(2)];
imshow(img2_rgb);
rectangle('Position',[window(1),window(2),col1,row1],'edgecolor','y');

%%
% More effient way: for example we choose the second frame here
window = [330,   213,   417,   261];
frame = 5;
fname = sprintf('%s\\frame00%d.jpg',path_to_images,frame+302);
img2_rgb = imread(fname);
img2_gray = rgb2gray(img2_rgb);
% smooth use gaussian kernel --sigma
gaussFilt = fspecial('gaussian',10,sigma);  % SEE 'help fspecial' 
img2_gray = imfilter(img2_gray, gaussFilt, 'symmetric'); 
img2_gray = im2double(img2_gray);

p = [0; 0];
step = 0;
threshold = 0.01;
max_iteration = 100;

[Tx, Ty] = meshgrid((1:col1)+window(1)-1, (1:row1)+window(2)-1);
[Gx, Gy] = imgradientxy(img1_gray, 'central');
Ix = Gx(:);
Iy = Gy(:);
A = [Ix, Iy];
for i=1:max_iteration
    % 1. Warp I with W
    warped = interp2(img2_gray, Tx + p(1), Ty + p(2), 'Bilinear');
    warped(isnan(warped)) = 0;
    % 2. Subtract I from T
    I_error = img1_gray(:) - warped(:);
    % 3. Multiply steepest descend with error and 8. compute dp
    dp =  pinv(A) * I_error;
    % 4. Update the parameters p <- p + delta_p
    p = p + dp;
    if norm(dp) < threshold
        break;
    end
    step = step + 1;
end
window = window + [p(1), p(2), p(1), p(2)];
imshow(img2_rgb);
rectangle('Position',[window(1),window(2),col1,row1],'edgecolor','y');