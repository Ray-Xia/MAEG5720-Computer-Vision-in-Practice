function [stitchedImage] = main_multiple_images(img1,img2,img3)

% clc;clear;
% img1 = imread('treeL.jpg');
% img2 = imread('treeC.jpg');
% img3 = imread('treeR.jpg');

img1 = imresize(img1, 0.5); img1 = rgb2gray(img1); img1 = single(img1);
img2 = imresize(img2, 0.5); img2 = rgb2gray(img2); img2 = single(img2); % reference image
img3 = imresize(img3, 0.5); img3 = rgb2gray(img3); img3 = single(img3);

[f1,d1] = vl_sift(img1); 
[f2,d2] = vl_sift(img2);
[f3,d3] = vl_sift(img3);

[matches12,scores12] = match_descriptor(d2,d1);
[bestH12,num_of_inliers12] = RANSAC(f2,f1,matches12);

[matches23,scores23] = match_descriptor(d2,d3);
[bestH23,num_of_inliers23] = RANSAC(f2,f3,matches23);

rows1 = size(img1,1); cols1 = size(img1,2);
rows2 = size(img2,1); cols2 = size(img2,2);
rows3 = size(img3,1); cols3 = size(img3,2);

stitchedImage = img2;
stitchedImage = padarray(stitchedImage, [0 cols1], 0,'pre');% 向左边填充
stitchedImage = padarray(stitchedImage, [0 cols3], 0,'post');% 向右边填充
stitchedImage = padarray(stitchedImage, [rows2 0], 0,'both');% 上下两边行填充
rows = size(stitchedImage,1);
cols = size(stitchedImage,2);

for i = 1:cols
    for j = 1:rows
            p1 = bestH12 * [i-cols1; j-rows1; 1];            
            p1 = p1 ./ p1(3);
            x1 = round(p1(1));
            y1 = round(p1(2));

            p3 = bestH23 * [i-cols3; j-rows3; 1];
            p3 = p3 ./ p3(3);
            x3 = round(p3(1));
            y3 = round(p3(2));
    
            if x1 > 0 && x1 <= cols1 && y1 > 0 && y1 <= rows1               
                    stitchedImage(j, i) = max(img1(y1, x1),stitchedImage(j, i));                                  
            end

            if x3 > 0 && x3 <= cols3 && y3 > 0 && y3 <= rows3               
                    stitchedImage(j, i) = max(img3(y3, x3),stitchedImage(j, i));                                  
            end
    end
end

% figure;imshow(stitchedImage,[])
stitchedImage = uint8(stitchedImage);

end