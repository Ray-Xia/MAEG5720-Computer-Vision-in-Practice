function [out] = prewitt_filter(img,threshold)
% Example: out = prewit_filter(img,threshold)
% Illustration: out is the output image; img is the input image; 
%               threshold is the threshold value defined by users

    input_image = rgb2gray(img);
    Mx = [-1 0 1; -1 0 1; -1 0 1]; 
    My = [-1 -1 -1; 0 0 0; 1 1 1];    
    
%     Gx = imfilter(input_image,Mx);
%     Gy = imfilter(input_image, My);
    Gx = box_filter(input_image,Mx);
    Gy = box_filter(input_image, My);
    filtered_image = uint8(sqrt(double(Gx.^2)+double(Gy.^2)));      

    out = max(filtered_image, threshold);
    out(out == round(threshold)) = 0;
    out = imbinarize(uint8(out));
end