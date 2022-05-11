function [window,step] = trackTemplate(img1,img2,window,sigma)

windowInit = [0,9,9,0];
step = 0;
p = [0; 0];    
threshold = 0.01;
max_iteration = 100;

img1_gray = rgb2gray(img1);
img2_gray = rgb2gray(img2);            
img1_gray = im2double(img1_gray);
img2_gray = im2double(img2_gray);
[row1, col1] = size(img1_gray);
[row2, col2] = size(img2_gray);

% smooth use gaussian kernel --sigma
gaussFilt = fspecial('gaussian',10,sigma);  
img1_gray = imfilter(img1_gray, gaussFilt, 'symmetric'); 
img2_gray = imfilter(img2_gray, gaussFilt, 'symmetric'); 

if isequal(window,windowInit)
    % %%%%%%%%%%%%%%%%%% use Nomarlized Cross-Correlation(NCC) to match the first frame with the template     
    % Initial bbox info using NCC
    I_NCC=normxcorr2(img1_gray,img2_gray);
    [ypeak,xpeak]=find(I_NCC==max(I_NCC(:))); % bottom right pt: y2, x2
    window = [xpeak-col1+1,ypeak-row1+1,xpeak,ypeak]; % x1,y1,x2,y2
    imshow(img2);
    rectangle('Position',[window(1),window(2),col1,row1],'edgecolor','y');
else
    %% %%%%%%%%%%%%%%%%%% use Lucas Kanade to deal with following frames        
    
    [Tx, Ty] = meshgrid((1:col1)+window(1)-1, (1:row1)+window(2)-1);
    for i=1:max_iteration
        % 1. Warp I with W
        warped = interp2(img2_gray, Tx + p(1), Ty + p(2), 'Bilinear');
        warped(isnan(warped)) = 0;
    
        % 2. Subtract I from T
        I_error = img1_gray(:) - warped(:);
    
        % 3. Compute Gradient
        % [Gx, Gy] = imgradientxy(warped, 'central');
        [Gx, Gy] = imgradientxy(img1_gray, 'central');
        Ix = Gx(:);
        Iy = Gy(:);
    
        % 4. Evaluate the jacobian
        W = [1, 0; 0, 1];
    
        % 5. Compute the steepest descent
        A = [Ix, Iy] * W;
    
        % 6. Compute inverse Hessian
        invH = pinv(A' * A);
    
        % 7. Multiply the steepest descent with error and 8. compute dp
        dp =  invH * A' * I_error;
    
        % 9. Update the parameters p <- p + delta_p
        p = p + dp;
        first_p = p;
        if norm(dp) < threshold
            break;
        end

        step = step + 1;
    end
    
    window = window + [p(1), p(2), p(1), p(2)];
    imshow(img2);
    rectangle('Position',[window(1),window(2),col1,row1],'edgecolor','y');
end

end