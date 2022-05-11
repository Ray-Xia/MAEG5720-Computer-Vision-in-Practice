function [cropedImage] = stitch(img1,img2,bestH)

%% stitch image using the best H
% clc;clear;
% img1 = imread('cu1.JPG');
% img2 = imread('cu2.JPG');
% img1 = imresize(img1, 0.5); img1 = rgb2gray(img1); img1 = single(img1);
% img2 = imresize(img2, 0.5); img2 = rgb2gray(img2); img2 = single(img2);

rows1 = size(img1,1); cols1 = size(img1,2);
rows2 = size(img2,1); cols2 = size(img2,2);
stitchedImage = img1;
stitchedImage = padarray(stitchedImage, [0 cols2], 0, 'post');% 向右边填充
% stitchedImage = padarray(stitchedImage, cols2, 0, 'post');
% figure;imshow(uint8(stitchedImage));
stitchedImage = padarray(stitchedImage, [rows2 0], 0, 'both');
% stitchedImage = padarray(stitchedImage, rows2 , 0, 'both');
% figure;imshow(uint8(stitchedImage))
rows = size(stitchedImage,1);
cols = size(stitchedImage,2);

for i = 1:cols
    for j = 1:rows
%         if rows2 < j && j < (rows2 + rows1)
            p2 = bestH * [i; j-rows2; 1];
            p2 = p2 ./ p2(3);
            x2 = round(p2(1));
            y2 = round(p2(2));
    
            if x2 > 0 && x2 <= cols && y2 > 0 && y2 <= rows
                if x2 < cols2 && y2 < rows2
                    stitchedImage(j, i) = max(img2(y2, x2),stitchedImage(j, i));
                end                  
            end
%         end
    end
end

% figure;imshow(stitchedImage,[])

%% Crop Image
sum_rows = sum(stitchedImage.').';
for k = 1:length(sum_rows)
    if sum_rows(k+1) > 0 && sum_rows(k) == 0
        ymin = k;        
    end

    if sum_rows(k+1) == 0 && sum_rows(k) > 0
        ymax = k+1;
        break;
    end
end

% sum_cols = sum(stitchedImage);
% for w = 1:length(sum_cols)
%     if sum_cols(w+1) == 0 && sum_cols(w) > 0
%         xmax = w+1;
%         break;
%     end
% 
% end

min_num_of_zeros = rows;
for j = 1: cols    
    num_of_zeros = 0;
    for i = 1:rows
        if stitchedImage(i,j) == 0
           num_of_zeros = num_of_zeros +1; 
        end
    end

    if num_of_zeros < min_num_of_zeros
        min_num_of_zeros = num_of_zeros;
        xmax = j; % found the column with the fewest number of zeros
    end
end

xmin = 0;
width = xmax -xmin;
height = ymax -ymin;

cropedImage = imcrop(stitchedImage,[xmin,ymin,width,height]);
cropedImage = uint8(cropedImage);
% figure;imshow(cropedImage,[])

end