clc; clear
%%
% HW1 reference solution
% related CV functions:
% Q1: imread/rgb2gray/imrotate/imtranslate/
%     affine2d/imwrap
% Q2: imresize: parameters (nearest, bilinear, bicubic)
% 
% Q3:
close all; clear; clc;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Questin 1

image=imread('lena.jpg');
figure; imshow(image); title('Q1a Color Image (10%) ');

greyimage=rgb2gray(image);
figure; imshow(greyimage); title('Q1b Grey Image (10%)');

% no matter they do it using the color and gray.
RotImage=imrotate(image, 30);
figure; imshow(RotImage); title('Q1c Rotated Image (10%)');

% some students may use the original image, as long as they translated.
TransImage=imtranslate(RotImage, [10, 10]);
figure; imshow(TransImage); ('Q1d Rotated and Translated Image (10%)');

scale = 1.0;
angle = 30*pi/180;
tx = 10;
ty = 10;
sc = scale*cos(angle);
ss = scale*sin(angle);

T = [sc -ss 0;
    ss sc 0;
    tx ty 1];

A= affine2d(T);

Image2=imwarp(image, A);
figure; imshow(Image2); title('Q1e Warped Image (10%)');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Question 2

halfimageNN = imresize(image, 0.5, 'nearest');
ExpImageNN = imresize(halfimageNN, 2, 'nearest');
figure; imshow(ExpImageNN); title('Q2a1 Nearest without smoothing (7%)');

ExpImageBL = imresize(halfimageNN, 2, 'bilinear');
figure; imshow(ExpImageBL); title('Q2a2 bilinear without smoothing (7%)');

ExpImageBC = imresize(halfimageNN, 2, 'bicubic');
figure; imshow(ExpImageBC); title('Q2a3 bicubic without smoothing (7%)');




halfimageBL = imresize(image, 0.5, 'bilinear');

ExpImageBL = imresize(halfimageBL, 2, 'nearest');
figure; imshow(ExpImageBL); title('Q2b1 Nearest with smoothing (7%)');

ExpImageBL = imresize(halfimageBL, 2, 'bilinear');
figure; imshow(ExpImageBL); title('Q2b2 bilinear with smoothing (7%)');

ExpImageBC = imresize(halfimageBL, 2, 'bicubic');
figure; imshow(ExpImageBC); title('Q2b3 bicubic with smoothing (7%)');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% For Q3 they can use whatever method to do. Just show the Nearest
% Neighbour result and Bilinear interpolation result without using 
% "reimsize" is good. Of course, the bilinear result should looks better.
% sample answers
% downsample
img_rsz_NN=my_resize(image,'NN',0.5);
img_rsz_BL=my_resize(image,'Bilinear',0.5);
figure
subplot(121),subimage(img_rsz_NN),title('Q3 resize image with self-defined NN')
subplot(122),subimage(img_rsz_BL),title('Q3 resize image with self-defined Bilinear')
% upsample
img_rsz_NN=my_resize(image,'NN',2);
img_rsz_BL=my_resize(image,'Bilinear',2);
figure
subplot(121),subimage(img_rsz_NN),title('Q3 resize image with self-defined NN')
subplot(122),subimage(img_rsz_BL),title('Q3 resize image with self-defined Bilinear')