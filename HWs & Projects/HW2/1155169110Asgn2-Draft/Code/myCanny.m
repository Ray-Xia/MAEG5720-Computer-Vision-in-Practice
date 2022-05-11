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

figure,imshow(uint8(output)); title('Non Maximum Suppression without Normalization');
mag = output;
% % Normalize the ouput
% mag = NormalizeMatrix(output);
% figure; imshow(mag);
% title('Non Maximum Suppression with Normalization');

%%
% (f)
% Perform double thresholding to find strong edges and weak edges
highThresholdRatio = 0.19;
lowThresholdRatio = 0.48;

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

