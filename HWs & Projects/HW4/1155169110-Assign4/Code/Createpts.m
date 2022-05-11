% Test Algorithm to generate projection of 3D points into 2 views
clc;clear;
close all;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% create two camera intrinsic
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% for camera one
f1=256;
hx1=128;
hy1=128;
K1=[f1, 0, hx1; 0, f1, hy1; 0, 0, 1];

% for camera two
f2=256;
hx2=128;
hy2=128;
K2=[f2, 0, hx2; 0, f2, hy2; 0, 0, 1];

% Create Camera extrinic parameters
% For Camera 1
angx = 15/180*pi ; angy = 45/180*pi; angz = 0;
tx = -35; ty = 0; tz = -100;
Rx = [ 1 0 0; 0 cos(angx) -sin(angx); 0 sin(angx) cos(angx)];
Ry = [ cos(angy) 0 sin(angy); 0 1 0; -sin(angy) 0 cos(angy)];
Rz = [ cos(angz) -sin(angz) 0; sin(angz) cos(angz) 0; 0 0 1];

R1=Rx*Ry*Rz;
M1 = [R1 [tx;ty;tz]; 0 0 0 1]; %Camera 1 matrix

% Create Camera extrinic parameters
% For Camera 2
angx = 15/180*pi; angy = -45/180*pi; angz = 0;
tx = 20; ty = 0; tz = -100;

Rx = [ 1 0 0; 0 cos(angx) -sin(angx); 0 sin(angx) cos(angx)];
Ry = [ cos(angy) 0 sin(angy); 0 1 0; -sin(angy) 0 cos(angy)];
Rz = [ cos(angz) -sin(angz) 0; sin(angz) cos(angz) 0; 0 0 1];

R2=Rx*Ry*Rz;
M2 = [R2 [tx;ty;tz]; 0 0 0 1]; %Camera 2 matrix


%Create 3D Points on Cube   %Construct image points x'i and x"i
Point_3D=[0,  0,   0,  2,  2,  2, 4, 4, 4,  0,  0,    0,  2,  2,  2, 4, 4, 4,    0,   0,     0,  2,  2,  2, 4, 4, 4; ...   
        0,  2,   4,   0,  2, 4,   0,  2, 4,  0,  2, 4,   0,  2, 4,   0,  2, 4,    0,  2,  4,   0,  2, 4,   0,  2, 4; ...
        0,  0,    0,   0,   0,   0,   0,   0,   0, 2,  2,  2,  2,  2,  2,  2,  2,  2,  4, 4,  4, 4, 4, 4, 4, 4, 4; ...
        1,  1,    1,   1,   1,   1,   1,   1,   1,  1,   1,   1,   1,   1,   1,   1,   1,   1,    1,   1,    1,   1,   1,   1,   1,   1,   1];     

% plot3(Point_3D(1,:),Point_3D(2,:),Point_3D(3,:), 'o');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Project the points into image 1
M1_Ext=M1(1:3,:);

Proj_Matrix_1=K1*M1_Ext;
p1 = Proj_Matrix_1*Point_3D;
p1(1,:) = p1(1,:)./p1(3,:);
p1(2,:) = p1(2,:)./p1(3,:);
p1(3,:) = p1(3,:)./p1(3,:);
 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Project the points into image 2
M2_Ext=M2(1:3,:);

Proj_Matrix_2=K2*M2_Ext;
p2 = Proj_Matrix_2*Point_3D;
p2(1,:) = p2(1,:)./p2(3,:);
p2(2,:) = p2(2,:)./p2(3,:);
p2(3,:) = p2(3,:)./p2(3,:);

% Add some noise to the image points
sigma = 0.2;
p1(1:2,:) = p1(1:2,:) + sigma*randn(2,size(p1, 2));

% Add some noise to the image points
sigma = 0.2;
p2(1:2,:) = p2(1:2,:) + sigma*randn(2,size(p2, 2));hold on

%% (a)
%Call FfromEightPnts to calculate the F-Matrix;
fMatrix = FfromEightPnts(p1, p2); 
% fMatrix = WithoutNormFfromEightPnts(p1, p2); 
F = fMatrix

% Calculate FM through MATLAB function 
p1_m = p1(1:2,:).';
p2_m = p2(1:2,:).';
[matlab_FM,inliers] = estimateFundamentalMatrix(p1_m,p2_m,"Method",'Norm8Point');

% error = abs(F - matlab_FM)

%% (b)
% Display image 1
figure;
I1 = zeros(512,512);
subplot(1,2,1);
imshow(I1);
title('Left Image');
hold on;
% plot(p1(1,:), p1(2,:), 'g*');
plot(p1(1,:), p1(2,:), 'w*');hold on

% (c) Draw the epipole
epipole1 = null(F);
plot(epipole1(1)/epipole1(3),epipole1(2)/epipole1(3),'go');hold on

% Draw epipolar line on Image I 
% imshow(I1); 
l1 =  F.' * p2;
% l1 =  matlab_FM.'* p2;
points = lineToBorderPoints(l1.',size(I1));
line( points(:,[1,3])', points(:,[2,4])' )  

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Display image 2 
I2 = zeros(512,512);
subplot(1,2,2);
imshow(I2);hold on;
title('Right Image');
% plot(p2(1,:), p2(2,:), 'r*');
plot(p2(1,:), p2(2,:), 'w*'); 

% (c) Draw the epipole 
epipole2 = null(F.');
plot(epipole2(1)/epipole2(3),epipole2(2)/epipole2(3),'go');hold on;

% Draw epipolar line on Image II 
% imshow(I1); 
l2 =  F * p1;
% l2 =  matlab_FM * p1;
points = lineToBorderPoints(l2.',size(I2));
line( points(:,[1,3])', points(:,[2,4])' )  

sgtitle('Condition the FM to rank 2 with Normalization')
% sgtitle('Not Condition the FM to rank 2 with Normalization')
% sgtitle('Without Normalization')
