function [bestH,num_of_inliers] = RANSAC(f1,f2,matches)

PixelError = 1;
num_of_matches = length(matches);
num_of_iteration = 100;
num_pts_ranSample = 4;
max_inliers = 0;

for i = 1:num_of_iteration
    num_of_inliers = 0;
    chosen_pairs = randsample(num_of_matches,num_pts_ranSample); % randomly choose 4 pairs
    pair1 = matches(:,chosen_pairs(1));
    pair2 = matches(:,chosen_pairs(2));
    pair3 = matches(:,chosen_pairs(3));
    pair4 = matches(:,chosen_pairs(4));
    x1 = f1(1,pair1(1,:)); y1 = f1(2,pair1(1,:));
    x2 = f1(1,pair2(1,:)); y2 = f1(2,pair2(1,:));
    x3 = f1(1,pair3(1,:)); y3 = f1(2,pair3(1,:));
    x4 = f1(1,pair4(1,:)); y4 = f1(2,pair4(1,:));
    xp1 = f2(1,pair1(2,:)); yp1 = f2(2,pair1(2,:));
    xp2 = f2(1,pair2(2,:)); yp2 = f2(2,pair2(2,:));
    xp3 = f2(1,pair3(2,:)); yp3 = f2(2,pair3(2,:));
    xp4 = f2(1,pair4(2,:)); yp4 = f2(2,pair4(2,:));
    H_hypothesis = FindHomography(x1,y1,x2,y2,x3,y3,x4,y4, xp1,yp1,xp2,yp2,xp3,yp3,xp4,yp4); % Calculate Homography with 4 randomly chosen points

    for j = 1:num_of_matches
        pair = matches(:,j);
        pt = H_hypothesis * [f1(1:2,pair(1));1];        
        x_tf = pt(1)/pt(3);y_tf = pt(2)/pt(3);
        x_img2 = f2(1,pair(2));y_img2 = f2(2,pair(2));
        distance = norm([x_tf - x_img2, y_tf - y_img2],2); % the distance between point in image2 and the transformed point originally in image1
        if distance < PixelError
            num_of_inliers = num_of_inliers +1;
        end
    end    

    if num_of_inliers > max_inliers
        bestH = H_hypothesis;
        max_inliers = num_of_inliers;    
    end
     
end

end