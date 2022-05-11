function [out] = sobel_filter(img,threshold)
% Example: out = prewit_sobel(img,threshold)
% Illustration: out is the output image; img is the input image; 
%               threshold is the threshold value defined by users

    input_image = rgb2gray(img);
    Mx = [-1 0 1; -2 0 2; -1 0 1]; 
    My = [-1 -2 -1; 0 0 0; 1 2 1];   
    
    Gx = box_filter(input_image,Mx);
    Gy = box_filter(input_image, My);
    filtered_image = uint8(sqrt(double(Gx.^2)+double(Gy.^2)));      

    out = max(filtered_image, threshold);
    out(out == round(threshold)) = 0;
    out = imbinarize(uint8(out));
end