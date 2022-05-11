clear all; close all; clc;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% create the camera intrinsic
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% for camera one
f1=512;
hx1=256;
hy1=256;
K1=[f1, 0, hx1; 0, f1, hy1; 0, 0, 1]


% Create Camera extrinic parameters
% For Camera 1
angx = 0; angy = 0; angz = 0;
cx = 0; cy = 0; cz = 100;
Rx = [ 1 0 0; 0 cos(angx) -sin(angx); 0 sin(angx) cos(angx)];
Ry = [ cos(angy) 0 sin(angy); 0 1 0; -sin(angy) 0 cos(angy)];
Rz = [ cos(angz) -sin(angz) 0; sin(angz) cos(angz) 0; 0 0 1];

R1=Rx*Ry*Rz;
M1 = [R1 -R1*[cx;cy;cz]; 0 0 0 1]


Point_3D=[-10,  -10,   -10,  0,  0,  0, 10, 10, 10,   -10,  -10, -10,  0,  0,  0, 10, 10, 10,    -10,   -10,     -10,  0,  0,  0, 10, 10, 10; ...   
          -10,  0,   10,   -10,  0, 10,   -10,  0, 10,  -10,  0, 10,   -10,  0, 10,   -10,  0, 10,    -10,  0,  10,   -10,  0, 10,   -10,  0, 10; ...
          -10,  -10,    -10,   -10,   -10,   -10,   -10,   -10,   -10, 0,  0,  0,  0,  0,  0,  0,  0,  0,  10, 10,  10, 10, 10, 10, 10, 10, 10; ...
          1,  1,    1,   1,   1,   1,   1,   1,   1,  1,   1,   1,   1,   1,   1,   1,   1,   1,    1,   1,    1,   1,   1,   1,   1,   1,   1];     


plot3(Point_3D(1,:),Point_3D(2,:),Point_3D(3,:), 'o');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Project the points into image 1
M1_Ext=M1(1:3,:);

Proj_Matrix_1=K1*M1_Ext;
p1 = Proj_Matrix_1*Point_3D;
p1(1,:) = p1(1,:)./p1(3,:);
p1(2,:) = p1(2,:)./p1(3,:);
p1(3,:) = p1(3,:)./p1(3,:);
 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Display image 1
  figure; ax=axes;

  I = zeros(512,512);
  imshow(I);
  hold on;
  plot(p1(1,:), p1(2,:), 'g*');
  
  % Add some noise to the image points
   sigma = 0;
   p1(1:2,:) = p1(1:2,:) + sigma*randn(2,size(p1, 2));
   plot(p1(1,:), p1(2,:), 'w*');
   hold off;
   
   DecomposeCamera(M1(1:3,:));
   
   
   
