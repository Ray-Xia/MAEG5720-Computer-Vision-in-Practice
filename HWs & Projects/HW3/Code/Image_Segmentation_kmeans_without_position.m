% Written by Pratik Jain
% Subscribe me on YouTube
% https://www.youtube.com/PratikJainTutorials


clc;
clear;
close all;

%Threshold on No. of itterations
n = 20;

%Treshold on the cost J
n1= 5;

% Read the image and its size
a1 = imread('Image path here');

b = size(a1);
N1 = b(1)*b(2);

k= input("Enter the value of k: ");

a1 = double(a1);
%making a vector of image

x1 = [];
x2 = [];
x3 = [];
for i=1:b(2)
    x1 = [x1;a1(:,i,1)];
end

for i=1:b(2)
    x2 = [x2;a1(:,i,2)];
end

for i=1:b(2)
    x3 = [x3;a1(:,i,3)];
end

x = [x1 x2 x3];

%Initialise the K-means
Ix1 = 255*rand(1,k);
Ix2 = 255*rand(1,k);
Ix3 = 255*rand(1,k);


mean1 = zeros(k,3);

for i= 1:k
    mean1(i,:) = [Ix1(i),Ix2(i),Ix3(i)];
end

%plot the pixels in 3d space
figure;
scatter3(x(:,1),x(:,2),x(:,3),'.','b');
xlabel('r');
ylabel('g');
zlabel('b');

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
        if (temp == zeros(1,3))
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

disp(J1(i1-1));

% plot the data with the current means
figure;
scatter3(x(:,1),x(:,2),x(:,3),'.','b');
xlabel('r');
ylabel('g');
zlabel('b');

hold on
scatter3(mean1(:,1),mean1(:,2),mean1(:,3),'filled','r');


% plot pixels allocated to the means in different colors
figure;
for i = 1:k
scatter3(x(find(Rnk(:,i)),1),x(find(Rnk(:,i)),2),x(find(Rnk(:,i)),3));
hold on
end
hold on 
scatter3(mean1(:,1),mean1(:,2),mean1(:,3),'filled','r');
xlabel('r');
ylabel('g');
zlabel('b');

% plotting only means
% figure;
% scatter3(mean1(:,1),mean1(:,2),mean1(:,3),'filled','r');
% xlabel('r');
% ylabel('g');
% zlabel('b');

% Show the image
j2=1;
p = 0;
Out = zeros(b(1),b(2),3);
for i2 = 1:b(2)*b(1)
    [temp1,Itemp] = max(Rnk(i2,:));
    Out(i2-p,j2,1) = mean1(Itemp,1);
    Out(i2-p,j2,2) = mean1(Itemp,2);
    Out(i2-p,j2,3) = mean1(Itemp,3);
    if i2-p == b(1)
        p = p+b(1);
        j2 = j2+1;
    end
end

% Showing the Results

figure;
imshow(uint8(Out))


end

% Getting the Output segmented Image
j2=1;
p = 0;
Out = zeros(b(1),b(2),3);
for i2 = 1:b(2)*b(1)
    [temp1,Itemp] = max(Rnk(i2,:));
    Out(i2-p,j2,1) = mean1(Itemp,1);
    Out(i2-p,j2,2) = mean1(Itemp,2);
    Out(i2-p,j2,3) = mean1(Itemp,3);
    if i2-p == b(1)
        p = p+b(1);
        j2 = j2+1;
    end
end

% Showing the Results

imshow(uint8(Out))
figure;
imshow(uint8(a1))

    
    
    