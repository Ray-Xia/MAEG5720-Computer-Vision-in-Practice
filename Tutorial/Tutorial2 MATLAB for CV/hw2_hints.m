close all; clear; clc;
%% ---------Problem(1)
img=imread('lena.jpg');
figure; imshow(img); title('Original Image');
[m, n, p] = size(img);

% Adding Noise
noise_img=imnoise(img,'gaussian');
figure, imshow(noise_img); title('noisy Image');

% function padarray: pad the image for filtering
help padarray

%define smoothing filter
average_filter = 1/25 * ones(5, 5);
% visualize filter
imshow(average_filter)
% bascic operation with filter
img_crop = noise_img(1:5, 1:5, 1);
multiply_elementwise = double(img_crop) .* average_filter;
sum_column = sum(multiply_elementwise);
sum_all = sum(sum_column);
final_result = uint8(sum_all);

%% 
%subtracted image
sub_image = double(noise_img)-double(img);
imshow(sub_image)

% 
% %sharpened image
sharp_img = double(img)+double(sub_image);

%
% Hints for Frequency Filter
% calculate the fft of image (gray) and do shifting 
img_gray = rgb2gray(img);
FT_img = fftshift(fft2(double(img_gray)));

% some operations here
D = fspecial('gaussian', size(img_gray), 40);
H = D/max(D(:));%D(:)将D中的所有元素重构成一个列向量

% visualize the filter
imshow(H)
imshow(H, []) % [] will do automatic adjustment for optimal contrast display

G = H.*FT_img; 

% Calculate the inverse FFT
out_img = real(ifft2(double(G))); 
imshow(uint8(out_img)) % remember to convert to uint8 type
imshow(img_gray)

%% ---------Problem(2)
% prewitt + sobel filter
Mx = [-1 0 1; -1 0 1; -1 0 1]; 
My = [-1 -1 -1; 0 0 0; 1 1 1]; 
% Gx = box_filter(input_image, Mx);
% Gy = box_filter(input_image, My);
% filtered_image = sqrt(Gx.^2+Gy.^2);
% maybe a threshold operation needs to be done to filtered_image

%% ---------Problem(3)
% canny
% Hints1 for NMS(non-maximum suppression): using if to determine which case for angles
% Hints2 for 3f: follow the description in assignment (recursive programming).