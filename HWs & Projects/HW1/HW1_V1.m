%% ---------Problem(1)
%% question(a)
image_color=imread("lena.jpg");
figure('Name','color image','NumberTitle','off')
imshow(image_color)

%% question(b)
image_gray=rgb2gray(image_color);
figure('Name','gray image','NumberTitle','off')
imshow(image_gray)

%% question(c)
rotate_image=imrotate(image_gray,30);
figure('Name','rotate image','NumberTitle','off')
imshow(rotate_image)


%% TEST
alpha=30/180*pi; % pure rotation (+: counter-clockwise) based affine function
tform = affine2d([cos(alpha) -sin(alpha) 0; sin(alpha) cos(alpha) 0; 0 0 1]);
img_affine = imwarp(image_color,tform);
figure(2)
imshow(img_affine)

%% question(d)
tranlate_image=imtranslate(rotate_image,[10,10]);
figure('Name','translate image','NumberTitle','off')
imshow(tranlate_image)
%% TEST
tx=10;
ty=10;
tform = affine2d([1 0 tx; 0 1 ty; 0 0 1]);
img_affine = imwarp(image_color,tform);
figure(3)
imshow(img_affine)

%% question(e)
alpha=30/180*pi;
Tx=10;
Ty=10;
T=[cos(alpha) -sin(alpha) 0; sin(alpha) cos(alpha) 0; Tx Ty 1];
tform=affine2d(T);
warp_image=imwarp(image_color,tform);
figure('Name','warp image','NumberTitle','off')
imshow(warp_image)
%%
alpha=30/180*pi;
Tx=100;
Ty=100;
R=[cos(alpha) -sin(alpha) 0; sin(alpha) cos(alpha) 0; 0 0 1];
T = [1 0 0;0 1 0;Tx Ty 1];
tform_rt = affine2d(R*T);
image=imwarp(image_color,tform_rt);
imshow(image)

%%
%平移变换
I = imread('lena.jpg');
tform = affine2d([1 0 0;0 1 0;200 200 1]);
J = imwarp(I,tform);
imshow(J)

%% -----------Problem(2)

%% question(a)
image_color=imread("lena.jpg");

image_down_nearest=imresize(image_color,0.5, 'nearest');
image_up_nearest=imresize(image_down_nearest,2, 'nearest');
image_up_bilinear=imresize(image_down_nearest,2, 'bilinear');
image_up_bicubic=imresize(image_down_nearest,2, 'bicubic');

figure('Name','Without filter, nearest','NumberTitle','off')
subplot(2,2,1), imshow(image_down_nearest),title('resolution×0.5, nearest')
subplot(2,2,2), imshow(image_up_nearest),title('resolution×2, nearest')
subplot(2,2,3), imshow(image_up_bilinear),title('resolution×2, bilinear')
subplot(2,2,4), imshow(image_up_bicubic),title('resolution×2, bicubic')

%% question(b)
image_color=imread("lena.jpg");

image_down_nearest=imresize(image_color,0.5, 'bilinear');
image_up_nearest=imresize(image_down_nearest,2, 'nearest');
image_up_bilinear=imresize(image_down_nearest,2, 'bilinear');
image_up_bicubic=imresize(image_down_nearest,2, 'bicubic');

figure('Name','With filter, bilinear','NumberTitle','off')
subplot(2,2,1), imshow(image_down_nearest),title('resolution×0.5, nearest')
subplot(2,2,2), imshow(image_up_nearest),title('resolution×2, nearest')
subplot(2,2,3), imshow(image_up_bilinear),title('resolution×2, bilinear')
subplot(2,2,4), imshow(image_up_bicubic),title('resolution×2, bicubic')
