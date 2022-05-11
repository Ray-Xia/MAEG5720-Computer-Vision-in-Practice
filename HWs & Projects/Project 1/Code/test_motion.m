function test_motion(path_to_images, numimages)

clc; clear;
currentFolder = pwd;
path_to_images = [currentFolder, '\CarSequence'];

pics = dir([path_to_images,'\frame*']); % list all frame pictures
numimages = length(pics); % Total number of images

fname = sprintf('%s\\car_template.jpg',path_to_images);
img1_rgb = imread(fname); %  Read Template

W = zeros(numimages,4);% Preserve tracked coordinates of each frame
steps_of_each_frame = zeros(10,numimages);
% first_p_of_each_frame = zeros(10,numimages);

%%%    To Do    %%%%%
% try different sigmas
for i =1:10
    sigma = i;
    %%%%%%%%%%%%%%%%%%%%%
    window = [0,9,9,0]; %x1 y1 x2 y2, the random Initial coordinates of top-left and bottom-right corners
    figure;
    for frame = 1:numimages
            % Reads next image in sequence
            fname = sprintf('%s\\frame00%d.jpg',path_to_images,frame+302);
            img2_rgb = imread(fname); % Read a frame
    
            %%%%%%%%%%%%%%%%%%%%%% To Do %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%     
    
            [new_window,step] = trackTemplate(img1_rgb,img2_rgb,window,sigma);
            steps_of_each_frame(i,frame) = step;
            str = 'The iteration number for frame %d is %d \n';
            fprintf(str,frame,step)
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
    %         current_frame = img2_rgb;
    %         imshow(current_frame);
    %         rectangle('Position',[new_window(3),new_window(1),88,49],'edgecolor','y');
    
            pause(0.01);
            W(frame,:) = new_window;
            window = new_window;
            
    end
    
    save 'track_coordinates' W;

end % end for

figure(111);
p = plot(1:numimages,steps_of_each_frame,'LineWidth',1);
p(1).Marker = 'o';
p(1).Color = 'g';
p(3).Marker = 'd';
p(3).Color = 'm';
p(5).Marker = '*';
p(5).Color = 'c';
p(7).Marker = 'x';
p(7).Color = 'r';
p(9).Marker = 's';
p(9).Color = 'k';
legend
xlim([1,110])
xlabel("Frame Number")
ylabel("Iteration Number")
title("Change in Interation Number of Each Frame When \sigma Changes")

end % end function
