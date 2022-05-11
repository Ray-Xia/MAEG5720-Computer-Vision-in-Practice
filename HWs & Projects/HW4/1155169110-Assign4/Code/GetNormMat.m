function Norm_Matrix=GetNormMat(pts)
% Function: Compute the normalized matrix of 2D points in HC
% Move the centroid to the origin
% Average distance of points from centroid is sqrt(2)

%convert the HC to Euclidean coordinates
pts(1,:) = pts(1,:) ./ pts(3,:); % x-coordinate
pts(2,:) = pts(2,:) ./ pts(3,:); % y-coordinate
pts = pts(1:2,:); % Now, each row is a Euclidean coordinate 

%calculate the centroid
ctrd = mean(pts,2);

%find the distance of each point from centroid
dist = sqrt(sum(pts-ctrd).^2);
%find the average distance to centroid
meanDist = mean(dist);

%scale it to sqrt(2);
s = sqrt(2) / meanDist;

%create the 3x3 matrix and return
Norm_Matrix = [s  0  -s*ctrd(1);
               0  s  -s*ctrd(2);
               0  0   1;];
 
% pts = [0 1 1 2; 2 5 2 3; 1 3 2 4]

end
