%----The following is the answer to Problem(1)
%% a)
image=imread('lena.jpg');
noisy_image=imnoise(image,'salt & pepper',0.1);
average_filter=ones(5,5)/25;%
noise_free=imfilter(noisy_image,average_filter);

subplot(2,2,1),imshow(image),title('Original Image');
subplot(2,2,2),imshow(noisy_image),title('Noisy Image');
subplot(2,2,3),imshow(noise_free),title('After Average Filtering');
%% b)
gray_image=rgb2gray(image);
b=box_filter(gray_image,average_filter);
figure;
imshow(uint8(b))
%% Test
clc;
clear;
close all;
xn=input('enter the sequence 1');
l1=length(xn);
hn=input('enter the sequence 2');
l2=length(hn);
m=l1+l2-1;

z_original=zeros(1,m);
z_original=conv(xn,hn);
disp(z_original);

z_new=zeros(1,m);

xn11=[xn,zeros(1,l2-1)];
disp(xn11);

hn11=[hn,zeros(1,l1-1)];
disp(hn11);

for i=1:m
    for j=1:i
        z_new(i)=z_new(i)+xn11(j)*hn11(i-j+1);
    end
end

disp(z_new);


%%
A = rand(3);
B = rand(4);
Cfull = conv2(A,B);
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
% imshow(uint8(img_gf1));

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

% subplot(2,2,4),imshow(uint8(img_diff3)),title('Difference 3');

% imshow(uint8(sharpend_img1));
% imshow(uint8(img_gf1));

%----The following is the answer to Problem(2)
%% (a)
% Please check the 'sobel_filter.m' for the sobel_filter function
% Here is a piece of test code for sobel_filter
clc;clear;
input_image = imread('lena.jpg');
out = prewitt_filter(input_image,100);
figure,imshow(out);title('Edge Detected Image with Prewitt Filter');
%%
input_image = imread('lena.jpg');
input_image = rgb2gray(input_image);
Mx = [-1 0 1; -1 0 1; -1 0 1]; 
My = [-1 -1 -1; 0 0 0; 1 1 1]; 
% img_prewitt=imfilter(input_image,Mx);
% figure,imshow(img_prewitt);

% Gx = imfilter(input_image,Mx);
% Gy = imfilter(input_image, My);
Gx = box_filter(input_image,Mx);
Gy = box_filter(input_image, My);
filtered_image = uint8(sqrt(double(Gx.^2)+double(Gy.^2)));
figure,imshow(filtered_image); title('Filtered Image');

% Define a threshold value
thresholdValue = 100; % varies between [0 255]
output_image = max(filtered_image, thresholdValue);
output_image(output_image == round(thresholdValue)) = 0;
figure,imshow(output_image); title('After Threshold');

% Displaying Output Image
output_image = imbinarize(uint8(output_image));
figure, imshow(output_image); title('Edge Detected Image');

%% (b)
% Please check the 'sobel_filter.m' for the sobel_filter function
% Here is a piece of test code for sobel_filter
clc;clear;
input_image = imread('lena.jpg');
out = sobel_filter(input_image,100);
figure,imshow(out);title('Edge Detected Image with Sobel Filter');
figure, imshow(output_image); title('Edge Detected Image');

%% ----The following is the answer to Problem(3)
close all;clear;clc;
% (a)
image=imread('lena.jpg');
image=rgb2gray(image);
figure, imshow(image); title('Gray Scale Image');


% (b)
gaussian_filter=fspecial('gaussian',[5,5],1);
image=imfilter(image,gaussian_filter);
figure, imshow(image); title('Image Filtered by Gaussian Filter');


% (c) Compute derivative using sobel filter
Mx = [-1 0 1; -2 0 2; -1 0 1]; 
My = [-1 -2 -1; 0 0 0; 1 2 1]; 
Gx = box_filter(image,Mx);
Gy = box_filter(image, My);


% (d) Compute magnitude and orientation
mag = sqrt(double(Gx.^2)+double(Gy.^2));
figure,imshow(uint8(mag)); title('Magnitude of the image');
% angle0 = atand(Gx./Gy);
angle = atan2d(Gy,Gx);
% angle1-angle0
% angle2 = atan2(Gy,Gx)*180/pi;
% angle2-angle1
% figure,imshow(uint8(angle)); title('Angle of the image');

%%
% (e) Non-maximum suppression
[h,w] = size(image);
% X=[-1,0,+1 ;-1,0,+1 ;-1,0,+1];
% Y=[-1,-1,-1 ;0,0,0 ;+1,+1,+1];    
output = zeros(h,w); % Initialize the output
% x = [0 1];
% 比较梯度方向的magnitude,找出极大值
for i=2:h-1 % row
        for j=2:w-1 % col         
            if (angle(i,j)>=-22.5 && angle(i,j)<=22.5) || ...
                (angle(i,j)<-157.5 && angle(i,j)>=-180)
                if (mag(i,j) >= mag(i,j+1)) && ...
                   (mag(i,j) >= mag(i,j-1))
                    output(i,j)= mag(i,j);
                else
                    output(i,j)=0;
                end
            elseif (angle(i,j)>=22.5 && angle(i,j)<=67.5) || ...
                (angle(i,j)<-112.5 && angle(i,j)>=-157.5)
                if (mag(i,j) >= mag(i+1,j+1)) && ...
                   (mag(i,j) >= mag(i-1,j-1))
                    output(i,j)= mag(i,j);
                else
                    output(i,j)=0;
                end
            elseif (angle(i,j)>=67.5 && angle(i,j)<=112.5) || ...
                (angle(i,j)<-67.5 && angle(i,j)>=-112.5)
                if (mag(i,j) >= mag(i+1,j)) && ...
                   (mag(i,j) >= mag(i-1,j))
                    output(i,j)= mag(i,j);
                else
                    output(i,j)=0;
                end
            elseif (angle(i,j)>=112.5 && angle(i,j)<=157.5) || ...
                (angle(i,j)<-22.5 && angle(i,j)>=-67.5)
                if (mag(i,j) >= mag(i+1,j-1)) && ...
                   (mag(i,j) >= mag(i-1,j+1))
                    output(i,j)= mag(i,j);
                else
                    output(i,j)=0;
                end
            end
        end
end

% figure,imshow(uint8(mag)); title('Non Maximum Suppression');
mag = NormalizeMatrix(output);
figure; imshow(mag);
title('Non Maximum Suppression');

%%
% (f)
% Perform double thresholding to find strong edges and weak edges
highThresholdRatio = 0.35;
lowThresholdRatio = 0.18;

highThreshold = max(max(mag))*highThresholdRatio;
lowThreshold = highThreshold*lowThresholdRatio;
strongEdgesRow = zeros(1,h); % For keeping track of the strong edge row index
strongEdgesCol = zeros(1,w); % For keeping track of the strong edge col index
weakEdgesRow = zeros(1,h);  % For keeping track of the weak edge row index
weakEdgesCol = zeros(1,w);  % For keeping track of the weak edge col index
strongIndex = 1;
weakIndex = 1;
for i=2:h-1 % row
    for j=2:w-1 % col
        if mag(i,j) > highThreshold    % Strong edge
            mag(i,j) = 1;
            strongEdgesRow(strongIndex) = i;
            strongEdgesCol(strongIndex) = j;
            strongIndex = strongIndex + 1;
        elseif mag(i,j) < lowThreshold % Not an edge
            mag(i,j) = 0;
        else                            % Weak edge
            weakEdgesRow(weakIndex) = i;
            weakEdgesCol(weakIndex) = j;
            weakIndex = weakIndex + 1;
        end
    end
end
figure; imshow(mag);
title('After Double Thresholding'); 

%%
% Determine actual edges in weak edges
set(0,'RecursionLimit',10000) % Set the maximum number of recursion
for i=1:strongIndex-1
    % Find the weak edges that are connected to strong edges and set 
    % them to 1
    mag = find_connected_weak_edge(mag, strongEdgesRow(i),...
        strongEdgesCol(i));
end
figure; imshow(mag);
title('Found Actual Edges from Weak Edges'); 

% Remove the remaining weak edges that are not actually edges
% and is noise instead
for i=1:weakIndex-1
    if mag(weakEdgesRow(i),weakEdgesCol(i)) ~= 1
        mag(weakEdgesRow(i),weakEdgesCol(i)) = 0;
    end
end
figure; imshow(mag);
title('Remaining Actual Edges'); 

%% --------------functions used-------------
% Normalize matrix
function[A] = NormalizeMatrix(A)
    A = A/max(A(:));
end

% Find weak edges that are connected to strong edges and set them to 1
function[mag] = find_connected_weak_edge(mag, r, c)
    for i = -2:1:2
        for j = -2:1:2
            if (r+i > 0) && (c+j > 0) && (r+i < size(mag,1)) && ...
                    (c+j < size(mag,2)) % Make sure we are not out of bounds
                if (mag(r+i,c+j) > 0) && (mag(r+i,c+j) < 1) % To find connected weak edge and set it to 1
                    mag(r+i,c+j) = 1;
                    mag = find_connected_weak_edge(mag, r+i, c+j);
                end
            end
        end
    end
end




