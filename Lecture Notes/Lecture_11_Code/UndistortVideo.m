% Pass in the camera parameters varibale from a saved work space
load MyUsbCamPara.mat 

% opens up user interface to select video file
filename = uigetfile('', 'select source video file');

% gets input file name
InputFileName = filename; 

% prompts the user to name the new video file
NewFileName = uiputfile; 

% Video file reader variable 
InputVideo = VideoReader(InputFileName);
OutputVideo = VideoWriter(NewFileName);
OutputVideo.FrameRate = InputVideo.FrameRate;
open(OutputVideo);


while hasFrame(InputVideo)
    Distorted_img = readFrame(InputVideo);
    UnDistorted = undistortImage(Distorted_img , cameraParams);
    writeVideo(OutputVideo, UnDistorted)
end

close(OutputVideo)
msgbox('conversion Complete')