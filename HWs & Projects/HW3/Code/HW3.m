clc;clear;close all;

k = 4; % number of clusters
threshold = 1; % Threshold for centroids change
% labels = 1:k;

% Read an RGB image and flatten it
rgbImage=imread('peppers.jpg');
subplot(1,2,1);
imshow(rgbImage);
numberOfPixels = size(rgbImage,1)*size(rgbImage,2);
% imflat = double(reshape(rgbImage,numberOfPixels,3));

%% Test for k_means()
[imgChange,labelChange] = k_means(rgbImage,10,15);
subplot(1,2,2);
imshow(uint8(imgChange));

%% Test for my_kmeans_function()
[clustered,C] = my_kmeans_function(imflat,4);
C = C/255;
clusteredImg = label2rgb(reshape(clustered,size(rgbImage,1),size(rgbImage,2)),C);
subplot(1,2,2);
imshow(clusteredImg);

%%
redChannel=rgbImage(:, :, 1);
greenChannel=rgbImage(:, :, 2);
blueChannel=rgbImage(:, :, 3);
imflat = double([redChannel(:), greenChannel(:), blueChannel(:)]); %将图像数据的每一个通道的所有像素点搞成一列
% % 知识点：A(:) 将 A 中的所有元素重构成一个列向量。
% % 如果 A 已经是列向量，则此表达式没有任何作用。

% Randomly select k points as initial centers
% R = randsample(imflat(:,1),k);
% G = randsample(imflat(:,2),k);
% B = randsample(imflat(:,3),k);

idxOfInitialCenters = randsample(numberOfPixels,k);% Index of Initial Centers
Centroids = imflat(idxOfInitialCenters,:);% Initial centroids
% Centroids = randi([1,255],k,3);
newCentroids = zeros(k,3);
labelOfEachPixel = zeros(numberOfPixels,1);% This vector contains label of each pixel
iterationNum = 0;
change = 1;

%%
% while iterationNum < 5
while change
    % Assign pixels to their nearest centroids.
    for i = 1:numberOfPixels  
%         minDist = inf;
        minDist = 100000;
        labelOfEachPixel(i) = 0;
        for j = 1:k
            difference = imflat(i,:)-Centroids(j);
            square = difference.^2;
            summation = sum(square);
            dist = sqrt(summation); % Get the distance from point i to centroid j
            if dist < minDist
                labelOfEachPixel(i) = j;
                minDist = dist;
            end
        end
    end

    % Calculate new means of each cluster   
    for x = 1:k
        totalNumOfx = 1;
        sumOfAllPixelValue = [0,0,0]; 
        for i = 1:numberOfPixels                           
            if labelOfEachPixel(i) == x
                sumOfAllPixelValue = sumOfAllPixelValue + imflat(i,:); 
                totalNumOfx = totalNumOfx + 1;
            end           

        end  
        newMean = sumOfAllPixelValue / totalNumOfx;
        newCentroids(x,:) = newMean;
%         idxOfLabelx = find(labelOfEachPixel==x); % Find all the index of Label x        
%         newCentroids(x,:) = mean(imflat(idxOfLabelx,:));% update centroids   
            
    end

    iterationNum = iterationNum +1; % Calculate the total iteration number

    % 结束循环条件
    diffCentrnoids = newCentroids -Centroids; % Difference between new Centroids and old Centroids
    normOfDiffCentrnoids = norm(diffCentrnoids);
%     comparisonResult = (diffCentrnoids < threshold); % Compare the difference with threshold
    if normOfDiffCentrnoids < threshold
%     if mean(reshape(comparisonResult,[],1)) == 1
        change = 0;
        % break;
    end
    Centroids = newCentroids; 

end

%%
% 图片显示
labelOfEachPixel = reshape(labelOfEachPixel,size(rgbImage,1),size(rgbImage,2));
Centroids = Centroids./255;
clusteredImage=label2rgb(labelOfEachPixel,Centroids);
subplot(1,2,2);
imshow(clusteredImage);
% reshape()

%%
% %迭代完了就更新图像
% for j = 1:k   %对所有簇心
%     for i = 1:numberOfPixels  %对所有像素点
%         if(labelOfEachPixel(i) == j)   %对应编号下的像素点
%            redChannel(i) = Centroids(j,1);
%            greenChannel(i) = Centroids(j,2);  %各个像素通道
%            blueChannel(i) = Centroids(j,3);
%         end     
%     end
% end
% 
% RGB_R = reshape(A,size_x,size_y);
% RGB_G = reshape(B,size_x,size_y);
% RGB_B = reshape(C,size_x,size_y);
% 
% image_new = pic;  %新图像
% image_new(:,:,1) = RGB_R;
% image_new(:,:,2) = RGB_G;
% image_new(:,:,3) = RGB_B;

%% kmeans() Test 1
% https://www.youtube.com/watch?v=okXd_ekO7pM

clc;clear;close all;
warning off

% Define K
numberOfClasses = 4;

rgbImage=imread('peppers.jpg');
subplot(1,2,1);
imshow(rgbImage);
redChannel=rgbImage(:, :, 1);
greenChannel=rgbImage(:, :, 2);
blueChannel=rgbImage(:, :, 3);
data=double([redChannel(:), greenChannel(:), blueChannel(:)]); %将图像数据的每一个通道的所有像素点搞成一列

[idx,C]=kmeans(data,numberOfClasses);% idx是每个观测值的簇索引的 n×1 向量 ，在 k×p 矩阵 C 中返回 k 个簇质心的位置。

idx=reshape(idx,size(rgbImage,1),size(rgbImage,2));
C=C/255;
clusteredImage=label2rgb(idx,C);
subplot(1,2,2);
imshow(clusteredImage);

%% kmeans() Test 2
% https://www.youtube.com/watch?v=okXd_ekO7pM

clc;clear;close all;
warning off
rgbImage = imread('DogNCat.jpg');
subplot(1,2,1);
imshow(rgbImage);
title('Original Image');
redChannel=rgbImage(:, :, 1);
greenChannel=rgbImage(:, :, 2);
blueChannel=rgbImage(:, :, 3);
data=double([redChannel(:), greenChannel(:), blueChannel(:)]);
for i=1:10
    numberOfClasses=i;
    [idx,C]=kmeans(data,numberOfClasses);
    idx=reshape(idx,size(rgbImage,1),size(rgbImage,2));
    C=C/255;
    clusteredImage=label2rgb(idx,C);
    subplot(1,2,2);
    imshow(clusteredImage);
    title(i);
    pause;
end
