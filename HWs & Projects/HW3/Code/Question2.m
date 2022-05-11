clc;clear;close all;

k = 4; % number of clusters
thCenChange = 0.25; % Threshold for centroids changed, range 0~1
thNumLabChange = 1;% Threshold for number of label changed

% Read an RGB image and flatten it
rgbImage=imread('lena.jpg');
[rows,cols,channels] = size(rgbImage);
subplot(1,2,1);
imshow(rgbImage);title('Original Image');
numberOfPixels = size(rgbImage,1)*size(rgbImage,2);
redChannel=rgbImage(:, :, 1);
redChannel=redChannel(:);
greenChannel=rgbImage(:, :, 2);
greenChannel=greenChannel(:);
blueChannel=rgbImage(:, :, 3);
blueChannel=blueChannel(:);
imflat = double([redChannel, greenChannel, blueChannel,zeros(numberOfPixels,1),zeros(numberOfPixels,1)]);
for j = 1:1:cols
    for i = 1:1:rows
        imflat(i+(j-1)*rows,4) = i;
        imflat(i+(j-1)*rows,5) = j;
    end
end

idxOfInitialCenters = randsample(numberOfPixels,k);% Index of Random Initial Centroids
Centroids = imflat(idxOfInitialCenters,:);% Initial centroids
oldCentroids = zeros(k,5);% Preserve the old centroids
labelOfEachPixel = zeros(numberOfPixels,1);% This vector contains label of each pixel
oldLabelOfEachPixel = zeros(numberOfPixels,1);% Preserve the old labels
iterationNum = 0; % Iteration number
change = 1;

while change
    % Assign pixels to their nearest centroids.
    for i = 1:numberOfPixels  
        minDist = 100000;
        labelOfEachPixel(i) = 0;
        for j = 1:k
            difference = imflat(i,:) - Centroids(j);
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
        totalNumOfx = 0;
        sumOfAllPixelValue = zeros(1,5); 
        for i = 1:numberOfPixels                           
            if labelOfEachPixel(i) == x
                sumOfAllPixelValue = sumOfAllPixelValue + imflat(i,:); 
                totalNumOfx = totalNumOfx + 1;
            end          
        end  
        fprintf('Label %d has %d in total.\n',x,totalNumOfx);
        newMean = sumOfAllPixelValue / totalNumOfx;
        Centroids(x,:) = newMean;            
    end

    iterationNum = iterationNum +1; % Calculate the total iteration number

    % End while for d(i)
    diffCentrnoids = Centroids - oldCentroids; % Difference between new Centroids and old Centroids
    normOfDiffCentrnoids = norm(diffCentrnoids);% Calculate the norm of diffCentrnoids
    ratio = normOfDiffCentrnoids / sqrt(255^2*3 + rows^2 + cols^2);
    if ratio < thCenChange
        change = 0;
    end

%     % End while for d(ii)
%     diffOfLabels = labelOfEachPixel - oldLabelOfEachPixel;
%     normOfDiffOfLabels = norm(diffOfLabels);
%     if normOfDiffOfLabels < thNumLabChange
%         change = 0;
%     end

    oldCentroids = Centroids; 
    oldLabelOfEachPixel = labelOfEachPixel;

end

% Show the clustered image
labelOfEachPixel = reshape(labelOfEachPixel,size(rgbImage,1),size(rgbImage,2));
Centroids = Centroids./255;
clusteredImage=label2rgb(labelOfEachPixel,Centroids(:,1:3));
subplot(1,2,2);
imshow(clusteredImage);title(['Clustered Image with K = ',num2str(k),', thCenChange = ',num2str(thCenChange)]);