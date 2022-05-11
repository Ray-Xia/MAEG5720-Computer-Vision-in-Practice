close all; clear; clc;
%----The following is the answer to Problem(1)
%% a)
image=imread('lena.jpg');
noisy_image=imnoise(image,'salt & pepper',0.1);
average_filter=ones(5,5)/25;
noise_free=imfilter(noisy_image,average_filter);

subplot(2,2,1),imshow(image),title('Original Image');
subplot(2,2,2),imshow(noisy_image),title('Noisy Image');
subplot(2,2,3),imshow(noise_free),title('After Average Filtering');
%% b)
% The following are testing codes. See the 'box_filter.m' file
img_avrg=box_filter(rgb2gray(image),average_filter);
figure,imshow(uint8(img_avrg)),title('Image filtered by average filter using my function');
%% c)
noisy_image_gray=rgb2gray(noisy_image);
gaussian_filter1=1/16*[1 2 1;2 4 2;1 2 1;];
gaussian_filter2=1/273*[1 4 7 4 1;4 16 26 16 4;7 26 41 26 7;4 16 26 16 4;1 4 7 4 1;];
gaussian_filter3=1/1003*[0 0 1 2 1 0 0;0 3 13 22 13 3 0;1 13 59 97 59 13 1;2 22 97 159 94 22 2;1 13 59 97 59 13 1;0 3 13 22 13 3 0;0 0 1 2 1 0 0;];

img_gf1=box_filter(noisy_image_gray,gaussian_filter1);
img_gf2=box_filter(noisy_image_gray,gaussian_filter2);
img_gf3=box_filter(noisy_image_gray,gaussian_filter3);

subplot(2,2,1),imshow(noisy_image_gray),title('Noisy Image');
subplot(2,2,2),imshow(uint8(img_gf1)),title('Noisy image filtered by gaussian filter1');
subplot(2,2,3),imshow(uint8(img_gf2)),title('Noisy image filtered by gaussian filter2');
subplot(2,2,4),imshow(uint8(img_gf3)),title('Noisy image filtered by gaussian filter3');

%% d)
img_diff1=double(noisy_image_gray)-double(img_gf1);
img_diff2=double(noisy_image_gray)-double(img_gf2);
img_diff3=double(noisy_image_gray)-double(img_gf3);

subplot(2,2,2),imshow(uint8(img_diff1)),title('Difference 1');
subplot(2,2,3),imshow(uint8(img_diff2)),title('Difference 2');
subplot(2,2,4),imshow(uint8(img_diff3)),title('Difference 3');

%% e)
sharpend_img1=double(img_diff1)+double(img_gf1);
sharpend_img2=double(img_diff2)+double(img_gf2);
sharpend_img3=double(img_diff1)+double(img_gf1);

subplot(2,2,1),imshow(uint8(sharpend_img1)),title('Sharpened image 1');
subplot(2,2,2),imshow(uint8(sharpend_img2)),title('Sharpened image 2');
subplot(2,2,3),imshow(uint8(sharpend_img3)),title('Sharpened image 3');
subplot(2,2,4),imshow(noisy_image_gray),title('noisy image gray');

%----The following is the answer to Problem(2)
%% (a)
% Please check the 'prewitt_filter.m' for the prewitt_filter function
% Here is a piece of test code for prewitt_filter
clc;clear;
input_image = imread('lena.jpg');
out = prewitt_filter(input_image,100);
figure,imshow(out);title('Edge Detected Image with Prewitt Filter');

%% (b)
% Please check the 'sobel_filter.m' for the sobel_filter function
% Here is a piece of test code for sobel_filter
clc;clear;
input_image = imread('lena.jpg');
out = sobel_filter(input_image,100);
figure,imshow(out);title('Edge Detected Image with Sobel Filter');

%%
%----The following is the answer to Problem(3)
% Please check the 'canny.m' for the canny function
% Here is a piece of test code for canny()
clear;close all;clc;
img=imread('lena.jpg');
canny(img,0.19,0.48,1);
