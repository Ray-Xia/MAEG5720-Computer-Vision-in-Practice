function test_motion(path_to_images, numimages)

clc; clear;

%%%    To Do    %%%%%
% try different sigmas
sigma = 10;
%%%%%%%%%%%%%%%%%%%%%
path_to_images = 'W:\CUHK\TA\2021R1_MAEG5720\Tutorials\Tutorial3 MATLAB for Tracking\CarSequence';
numimages = 101;
fname = sprintf('%s\\car_template.jpg',path_to_images);
img1_rgb = imread(fname);
img1_gray = rgb2gray(img1_rgb);
window = [0,9,9,0]; %x1 y1 x2 y2, the coordinates of top-left and bottom-right corners
W = zeros(numimages,4);
figure;
for frame = 1:numimages
        % Reads next image in sequence
        fname = sprintf('%s\\frame00%d.jpg',path_to_images,frame+302);
        img2_rgb = imread(fname);
        img2_gray = rgb2gray(img2_rgb);
        %%%%%%%%%%%%%%%%%%%%%% To Do %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        new_window = trackTemplate(img1_gray,img2_gray,window,sigma);
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        current_frame = img2_rgb;
        imshow(current_frame);
        rectangle('Position',[new_window(3),new_window(1),88,49],'edgecolor','y');

        pause(0.01);
        W(frame,:) = new_window;
        window = new_window;
        frame;
        
end
    save 'track_coordinates' W;
end
