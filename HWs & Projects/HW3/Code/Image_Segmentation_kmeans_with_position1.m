% Written by Pratik Jain
% Subscribe me on YouTube
% https://www.youtube.com/PratikJainTutorials

clc;
clear;
close all;
tic
%Threshold on No. of itterations
n = 20;

%Treshold on the cost J
n1= 0.1;

% Read the image and its size
a = imread('DogNCat.jpg');
a = double(a);
b = size(a);
N1 = b(1)*b(2);%first x second dimension of the image

k= input("Enter the value of k: ");%Input the number of clusters

%make vectors of image
x1 = [];
x2 = [];
x3 = [];
for i=1:b(2)
    x1 = [x1;a(:,i,1)];% link each column of first channel in series
end

for i=1:b(2)
    x2 = [x2;a(:,i,2)];% link each column of second channel in series
end

for i=1:b(2)
    x3 = [x3;a(:,i,3)];% link each column of third channel in series
end


x4 = zeros(1,b(1)*b(2));
x5 = zeros(1,b(1)*b(2));
p=1;
p1=0;
for v=1:b(2)
    for v1=1:b(1)
        x4(v1+p1)= v1;
        x5(v1+p1)= p;
    end
    p=p+1;
    p1 = p1+b(1);
end

x = [x1 x2 x3 x4' x5'];

%Initialise the K-means
Ix1 = 255*rand(1,k);
Ix2 = 255*rand(1,k);
Ix3 = 255*rand(1,k);
Ix4 = size(a,1)*rand(1,k);
Ix5 = size(a,2)*rand(1,k);


mean1 = zeros(k,5);

for i= 1:k
    mean1(i,:) = [Ix1(i),Ix2(i),Ix3(i),Ix4(i),Ix5(i)];
end

% Loop for alloting the samples to a mean and then reupdating the Means
i1=2;
J1 = [0 0];

while i1<n
    Rnk = zeros(N1,k);
    d = zeros(N1,k);
    for i=1:N1
        for j=1:k
            d(i,j) = (x(i,:)-mean1(j,:))*(x(i,:)-mean1(j,:))';
        end
        
        [min, Imin] = max(-d(i,:));
        Rnk(i,Imin) = 1;
        
        if sum(Rnk(i,:)) == 0
            po(i) = 1;
        end
    end
    
    J1(i1) = 0;
    sumRnk = zeros(1,k);
    for i=1:N1
        for j = 1:k
            J1(i1) = J1(i1)+Rnk(i,j)*d(i,j);
        end
        sumRnk = sumRnk + Rnk(i,:);
    end
    
    for i=1:N1
        for j = 1:k
            temp = Rnk(i,j)*x(i,:);
            if (temp == zeros(1,5))
            else
            mean1(j,:) = mean1(j,:)+temp;
            end
        end
    end
    for i =1:k
        if sumRnk(i)~=0
        mean1(i,:) = mean1(i,:)/sumRnk(i);
        end
    end
    
    
    if (abs(J1(i1)-J1(i1-1))<n1)
        break;
    end
    i1 = i1+1;
    J1 = [J1 0];
    
    disp(J1(i1-1));
    
    % Getting the Output segmented Image    
    j2=1;
    p = 0;
    Out = zeros(b(1),b(2));
    for i2 = 1:b(2)*b(1)
        [temp1,Itemp] = max(Rnk(i2,:));
        Out(i2-p,j2,1) = (mean1(Itemp,1));
        Out(i2-p,j2,2) = (mean1(Itemp,2));
        Out(i2-p,j2,3) = (mean1(Itemp,3));
        if i2-p == b(1)
            p = p+b(1);
            j2 = j2+1;
        end
    end    
    
end

figure;
imshow(uint8(Out));
%%
% Getting the Output segmented Image
color1 = linspace(0,255,k);
color = [color1' color1' color1'];
j2=1;
p = 0;
Out = zeros(b(1),b(2));
for i2 = 1:b(2)*b(1)
    [temp1,Itemp] = max(Rnk(i2,:));
    Out(i2-p,j2,1) = 1/3*(mean1(Itemp,1));
    Out(i2-p,j2,2) = 1/3*(mean1(Itemp,2));
    Out(i2-p,j2,3) = 1/3*(mean1(Itemp,3));

    if i2-p == b(1)
        p = p+b(1);
        j2 = j2+1;
    end
end

% Showing the Results

imshow(uint8(Out))

figure;
imshow(uint8(a))
toc
%%  
    
    