clc; clear
%%
%-------------------------------------------- using self-defined functions
% useful in large projects (simplify the code structure)
% below is the customized add function. (extra my_add.m file in the same folder)
c = my_add(1,2);

%% matrix operation (linear algebra)
%-------------------------------------------- build a mxn matrix (here we build a 3x3 matrix)
A = [1,-2,3; 
    4,-5,6;
    7,8,-9]; % 3x3 matrix build in function: eye(3), diag([2,3,4])
% How to take a number in a matrix? --using row and column index~
% element manner
A
A(1,3) % 1: the 1st row; 3: the 3th column
% column manner
A(1:2, 2) % 1:2 represent the 1st to 2nd row; 2: the 2nd column
% row manner
A(1, 1:3) % 1: the 1st row; 1:3 represent the 1st to 3th column
A(1, 1:end) % the same as the A(1, :) -- 1: the 1st row; 1:end represent the 1st to end column
det(A) % determinant
rank(A) % rank
inv(A) % inverse
pinv(A) % pseudo inverse
trace(A) % trace
eig(A) % eigen value
% be familar with using help function
help eig
E = eig(A)
[V,D] = eig(A)
% check the eigen values and corresponding eigen vectors (using subcoordinate to access the elements in a matrix)
% we can find the same results A*x = lamda * x
A*V(:, 1)
D(1,1)*V(:, 1)
% singular value decomposition
[U,S,V] = svd(A) % orthogonal
A
U*S*V'

% matrix built operation (NOTE: Avoid for loop in MATLAB, take advantage of parallel computing in MATLAB)
matA = randn(2,3); % matrix A with shape 2x3
matB = randn(3,4); % matrix B with shape 3x4
matC = randn(2,3); % matrix C with shape 2x3

matA * matB % matrix multiply operation --> return 2x4 matrix
matA .* matC % element mutiply --> return 2x3 matrix
matA .^ matC

% quick check
B = [2, 3; 4, 8; 2, 9];
C = [1,4; 2, 4; 1, 6];
B .^ C
B ./ C

%--------------------------------------------basic plot skills
% 2d plot
X = linspace(0, 2*pi, 100); % Equally spaced 100 points X within [0, 2*pi]
Y = sin(X); % parallel compute Y
plot(X, Y); % plot two variable
plot(Y); % plot one variable (Explore more options in PLOTS panel)
% 3d plot
[X,Y] = meshgrid(-2:.2:2, -2:.2:2);                                
Z = X .* exp(-X.^2 - Y.^2);                                        
surf(X,Y,Z)

%%
%-------------------------------------------- start with image processing
img_color = imread('test.jpg');
img_gray = rgb2gray(img_color);
figure(1)
subplot(1,2,1), subimage(img_color)
subplot(1,2,2), subimage(img_gray)

% crop a image with matrix operation
figure(2)
subplot(1,2,1), subimage(img_gray)
subplot(1,2,2), subimage(img_gray(1:128, 1:200)) % select a part of image to show

% basic image processing
% resize image
img_resize = imresize(img_color, 0.5);
figure(3)
subplot(1,2,1), subimage(img_color) % imshow(img_color)
subplot(1,2,2), subimage(img_resize) % imshow(img_resize)

img_resize = imresize(img_color, [500, 500]);
figure(4)
subplot(1,2,1), subimage(img_color) % imshow(img_color)
subplot(1,2,2), subimage(img_resize) % imshow(img_resize)

% flip image horizontally and vertically
img_flip_hz = img_color(:,end:-1:1,:);
figure(5)
subplot(1,2,1), subimage(img_color), title('original') % imshow
subplot(1,2,2), subimage(img_flip_hz), title('flip horizontally') % imshow

img_flip_vt = img_color(end:-1:1,:,:);
figure(6)
subplot(1,2,1), subimage(img_color), title('original') % imshow
subplot(1,2,2), subimage(img_flip_vt), title('flip vertically') % imshow

% affine
tform = affine2d([2 0.33 0; 0 1 0; 0 0 1]); % shear image
img_affine = imwarp(img_color,tform);
figure(7)
subplot(1,2,1), subimage(img_color), title('original') % imshow
subplot(1,2,2), subimage(img_affine), title('affined') % imshow

alpha=10/180*pi; % pure rotation (+: counter-clockwise) based affine function
tform = affine2d([cos(alpha) -sin(alpha) 0; sin(alpha) cos(alpha) 0; 0 0 1]);
img_affine = imwarp(img_color,tform);
figure(8)
subplot(1,2,1), subimage(img_color), title('original') % imshow
subplot(1,2,2), subimage(img_affine), title('affined') % imshow
%%
%-------------------------------------------- understanding HSV (hue, saturation, value(brightness))
% read more on http://color.lukas-stratmann.com/color-systems/hsv.html
ball = imread('colorball.jpg');
ball_hsv = rgb2hsv(ball);
ball_rgb = hsv2rgb(ball_hsv);
figure(9)
subplot(1,3,1), subimage(ball), title('original')
subplot(1,3,2), subimage(ball_hsv), title('hsv shown in rgb')
subplot(1,3,3), subimage(ball_rgb), title('rgb->hsv->rgb')

ball_hsv(:,:,1)=0.5*ball_hsv(:,:,1); % hue
% ball_hsv(:,:,2)=2.0*ball_hsv(:,:,2); % saturation
% ball_hsv(:,:,3)=0.5*ball_hsv(:,:,3); % value(brightness)
ball_ = hsv2rgb(ball_hsv);
figure(10)
subplot(1,2,1), subimage(ball), title('original')
subplot(1,2,2), subimage(ball_), title('operation in hsv')

%%
%-------------------------------------------- basic color segmentation demo using HSV
ball = imread('colorball.jpg');
ball_hsv = rgb2hsv(ball); % transfer the rgb channal into hsv channal
ball_new = 255*ones(size(ball)); % create a white (255 in all channals) template
ball_new_hsv = rgb2hsv(ball_new);
%------------ select the purple color
% [row, col] = ind2sub(size(ball_hsv), find(ball_hsv(:,:,1) > 125/180 ...
%     & ball_hsv(:,:,1) < 155/180 & ball_hsv(:,:,2) > 43/255 & ball_hsv(:,:,3) > 46/255));
%------------ select the yellow color
[row, col] = ind2sub(size(ball_hsv), find(ball_hsv(:,:,1) > 20/180 ...
    & ball_hsv(:,:,1) < 34/180 & ball_hsv(:,:,2) > 43/255 & ball_hsv(:,:,3) > 46/255));
% %------------ select the blue color
% [row, col] = ind2sub(size(ball_hsv), find(ball_hsv(:,:,1) > 35/180 ...
%     & ball_hsv(:,:,1) < 100/180 & ball_hsv(:,:,2) > 43/255 & ball_hsv(:,:,3) > 46/255));
for i = 1:length(row)
    ball_new_hsv(row(i), col(i), :) = ball_hsv(row(i), col(i), :);
end
ball_selected_color = hsv2rgb(ball_new_hsv); % To view the segmenatation results, we need transfer the hsv channal back to rgb channal
figure(11)
subplot(1,2,1), subimage(ball)
subplot(1,2,2), subimage(ball_selected_color)

%% computer vision toolbox (APPS->IMAGE PROCESSING AND COMPUTER VISION->Camera Calibrator)
% Demonstrate how to use calibration toolbox to calibrate cameras (experience: 15+ checkboard images with various position are needed)
% intrinsics = cameraParams.IntrinsicMatrix
% distortion = cameraParams.RadialDistortion