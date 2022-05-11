% Perform edge detection without interpolation during non maximum 
% suppression
function CannyEdgeDetectorNoInterpolation()
    close all;  % Close figures
    sigma = 1.4; % Gaussian filter sigma
    highThresholdRatio = 0.2; % High threshold ratio
    lowThresholdRatio = 0.15; % Low threshold ratio
    
    % Change the current folder to the folder of this m-file.
    % Courtesy of Brett Shoelson
    if(~isdeployed)
      cd(fileparts(which(mfilename)));
    end
    
    im = imread('lena.jpg');
    figure; imshow(im);
    title('Original Image');
    
    % Smooth with Gaussian 5x5 filter to reduce noise
    im = rgb2gray(im);
    figure; imshow(im);
    title('B/W Image');
    
    im = double(imgaussfilt(im,sigma));
    figure; imshow(NormalizeMatrix(im));
    title('Gaussian Filter');
    
    % Find the intensity gradient of the image
    Gx = SobelFilter(im, 'x');
    Gy = SobelFilter(im, 'y');
    Gx = imgaussfilt(Gx,sigma);
    Gy = imgaussfilt(Gy,sigma);
    figure; imshow(NormalizeMatrix(Gx));
    title('Gx Sobel Filter');
    figure; imshow(NormalizeMatrix(Gy));
    title('Gy Sobel Filter');
    
    % Find the magnitude of the gradient
    Gmag = sqrt(Gx.^2 + Gy.^2);
    angle = atan2(Gy,Gx)*180/pi;
    figure; imshow(NormalizeMatrix(Gmag));
    title('Gmag');
         
    % Perform non-maximum suppression without interpolation
    [h,w] = size(im);
    X=[-1,0,+1 ;-1,0,+1 ;-1,0,+1];
	Y=[-1,-1,-1 ;0,0,0 ;+1,+1,+1];
    output = zeros(h,w);
    x = [0 1];
    for i=2:h-1 % row
        for j=2:w-1 % col         
            if (angle(i,j)>=-22.5 && angle(i,j)<=22.5) || ...
                (angle(i,j)<-157.5 && angle(i,j)>=-180)
                if (Gmag(i,j) >= Gmag(i,j+1)) && ...
                   (Gmag(i,j) >= Gmag(i,j-1))
                    output(i,j)= Gmag(i,j);
                else
                    output(i,j)=0;
                end
            elseif (angle(i,j)>=22.5 && angle(i,j)<=67.5) || ...
                (angle(i,j)<-112.5 && angle(i,j)>=-157.5)
                if (Gmag(i,j) >= Gmag(i+1,j+1)) && ...
                   (Gmag(i,j) >= Gmag(i-1,j-1))
                    output(i,j)= Gmag(i,j);
                else
                    output(i,j)=0;
                end
            elseif (angle(i,j)>=67.5 && angle(i,j)<=112.5) || ...
                (angle(i,j)<-67.5 && angle(i,j)>=-112.5)
                if (Gmag(i,j) >= Gmag(i+1,j)) && ...
                   (Gmag(i,j) >= Gmag(i-1,j))
                    output(i,j)= Gmag(i,j);
                else
                    output(i,j)=0;
                end
            elseif (angle(i,j)>=112.5 && angle(i,j)<=157.5) || ...
                (angle(i,j)<-22.5 && angle(i,j)>=-67.5)
                if (Gmag(i,j) >= Gmag(i+1,j-1)) && ...
                   (Gmag(i,j) >= Gmag(i-1,j+1))
                    output(i,j)= Gmag(i,j);
                else
                    output(i,j)=0;
                end
            end
        end
    end
    
    Gmag = NormalizeMatrix(output);
    figure; imshow(Gmag);
    title('Non Maximum Suppression');
    
    % Perform double thresholding
    highThreshold = max(max(Gmag))*highThresholdRatio;
    lowThreshold = highThreshold*lowThresholdRatio;
    strongEdgesRow = zeros(1,h); % Keep track of the strong edge row index
    strongEdgesCol = zeros(1,w); % Keep track of the strong edge col index
    weakEdgesRow = zeros(1,h);  % Keep track of the weak edge row index
    weakEdgesCol = zeros(1,w);  % Keep track of the weak edge col index
    strongIndex = 1;
    weakIndex = 1;
    for i=2:h-1 % row
        for j=2:w-1 % col
            if Gmag(i,j) > highThreshold    % Strong edge
                Gmag(i,j) = 1;
                strongEdgesRow(strongIndex) = i;
                strongEdgesCol(strongIndex) = j;
                strongIndex = strongIndex + 1;
            elseif Gmag(i,j) < lowThreshold % No edge
                Gmag(i,j) = 0;
            else                            % Weak edge
                weakEdgesRow(weakIndex) = i;
                weakEdgesCol(weakIndex) = j;
                weakIndex = weakIndex + 1;
            end
        end
    end
    figure; imshow(Gmag);
    title('Double Threshold'); 
    
    % Perform edge tracking by hysteresis
    set(0,'RecursionLimit',10000)
    for i=1:strongIndex-1
        % Find the weak edges that are connected to strong edges and set 
        % them to 1
        Gmag = FindConnectedWeakEdges(Gmag, strongEdgesRow(i),...
            strongEdgesCol(i));
    end
    figure; imshow(Gmag);
    title('Edge Tracking Before Clean Up'); 
    
    % Remove the remaining weak edges that are not actually edges
    % and is noise instead
    for i=1:weakIndex-1
        if Gmag(weakEdgesRow(i),weakEdgesCol(i)) ~= 1
            Gmag(weakEdgesRow(i),weakEdgesCol(i)) = 0;
        end
    end
    figure; imshow(Gmag);
    title('Edge Tracking After Clean Up'); 
    
    % MATLAB canny comparison
    im = imread('lena.jpg');
    im = rgb2gray(im);
    im = edge(im, 'canny');
    figure; imshow(im);
    title('MATLAB');
end

% Normalize matrix
function[A] = NormalizeMatrix(A)
    A = A/max(A(:));
end

% Perform sobel filter
function[A] = SobelFilter(A, filterDirection)
    switch filterDirection
        case 'x' 
            Gx = [-1 0 +1; -2 0 +2; -1 0 +1];
            A = imfilter(A, double(Gx), 'conv', 'replicate');
        case 'y'
            Gy = [-1 -2 -1; 0 0 0; +1 +2 +1];
            A = imfilter(A, double(Gy), 'conv', 'replicate');
        otherwise
            error('Bad filter direction - try inputs ''x'' or ''y''');
    end
end

% Find weak edges that are connected to strong edges and set them to 1
function[Gmag] = FindConnectedWeakEdges(Gmag, row, col)
    for i = -2:1:2
        for j = -2:1:2
            if (row+i > 0) && (col+j > 0) && (row+i < size(Gmag,1)) && ...
                    (col+j < size(Gmag,2)) % Make sure we are not out of bounds
                if (Gmag(row+i,col+j) > 0) && (Gmag(row+i,col+j) < 1)
                    Gmag(row+i,col+j) = 1;
                    Gmag = FindConnectedWeakEdges(Gmag, row+i, col+j);
                end
            end
        end
    end
end