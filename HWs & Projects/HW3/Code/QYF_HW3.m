% Codes written by QIU Yufu
clc;clear;
pic = imread('peppers.jpg');   %读入图像

[size_x,size_y,dimension] = size(pic);   %得到图像大小

%定义簇心
K = 12;  %定义簇的个数

%-------------problem 1--------------%
% clusters = zeros(K, dimension); %定义簇心
% clusters(:,1:3) = randi([1,255],K,3);
%-------------problem 1--------------%

%-------------problem 2--------------%
clusters = zeros(K, dimension+2); %定义簇心
clusters(:,1:3) = randi([1,255],K,3);
clusters(:,4) = randperm(size_x,K);
clusters(:,5) = randperm(size_y,K);
%-------------problem 2--------------%

% 将图像RGB3通道分解
A = double(reshape(pic(:, :, 1), size_x * size_y, 1));   
B = double(reshape(pic(:, :, 2), size_x * size_y, 1));
C = double(reshape(pic(:, :, 3), size_x * size_y, 1));
xx = 1:1:size_x;
X = repmat((xx), 1, size_y);
yy = size_x * size_y;
Y = ones(1, yy);
for i = 1:yy
    Y(i) = Y(i) + floor(i/size_x);
end

sample_num = size(A,1);   %像素点个数
dist = zeros(1,K);

c = zeros(sample_num,1);
dist_min = 100000 ;
old_clusters = 1000;  %设置初值，上一次迭代的簇
old_c = 1000 * ones(sample_num,1);     %设置初值，上一次迭代的序列

%-------------problem 2--------------%
% kx = 0;   %x，y在分类中的权重，太大会导致分不出
% ky = 0;
%-------------problem 2--------------%

iter = 30; % 假定最多迭代30次
t = 0;
while(1)
%-------------c-------------%
%     t = t+1;  %(c)
%     if (t == iter)  %阈值
%        break;
%     end
%---------------------------%

%------------d(i)--------------%   
%     bb = clusters - old_clusters;
%     bbb = norm(bb);
%     t = t + 1;
%     if (bbb < 1)  %阈值
%         break;
%     end
%     old_clusters = clusters;
 %-----------------------------%   
 
 %------------d(ii)--------------%   
    cc = c - old_c;    %序列变化
    ccc = norm(cc);
    t = t + 1;
    if (ccc < 0.5)  %阈值
        break;
    end
    old_c = c;
 %-----------------------------%   
 
    for i = 1:sample_num  % 对于每个像素点
        dist_min = 100000;
        for j = 1:K   %对于每个簇cj
            % 计算向量模长
       %--------------------problem 1----------------------%
%             dist(j) = ((A(i) - clusters(j,1))^2 + (B(i) - clusters(j,2))^2 + (C(i) - clusters(j,3))^2)^0.5;  
       %--------------------problem 1----------------------%
          
       %--------------------problem 2----------------------%
            dist(j) = ((A(i) - clusters(j,1))^2 + (B(i) - clusters(j,2))^2 + (C(i) - clusters(j,3))^2 + ((X(i) - clusters(j,4))/size_x * 255).^2 + ((Y(i) - clusters(j,5))./size_y * 255).^2).^0.5;
       %--------------------problem 2----------------------%
            if(dist(j) < dist_min)
                c(i) = j;
                dist_min = dist(j);
                
            end 
        end
    end
    
    %更新簇心 
    for j = 1:K   %对所有簇心
        %初始化
        num = 0;
        A_sum = 0;
        B_sum = 0;
        C_sum = 0;
%------------problem 2-------------%
        X_sum = 0;
        Y_sum = 0;
%------------problem 2-------------%
        for i = 1:sample_num  %对所有像素点
            if(c(i) == j)  %找到对应编号下的像素点
                %求像素点均值准备
                num = num + 1;
                A_sum = A_sum + A(i);
                B_sum = B_sum + B(i);
                C_sum = C_sum + C(i);
            %------------problem 2-------------%
                X_sum = X_sum + X(i);
                Y_sum = Y_sum + Y(i);
            %------------problem 2-------------%
            end
        end
        %求像素点均值，赋值给簇心
        clusters(j,1) = A_sum/num;
        clusters(j,2) = B_sum/num;
        clusters(j,3) = C_sum/num;
    %------------problem 2-------------%
        clusters(j,4) = X_sum/num;
        clusters(j,5) = Y_sum/num;
    %------------problem 2-------------%
    end
end

%迭代完了就更新图像
for j = 1:K   %对所有簇心
    for i = 1:sample_num  %对所有像素点
        if(c(i) == j)   %对应编号下的像素点
           A(i) = clusters(j,1);
           B(i) = clusters(j,2);  %各个像素通道
           C(i) = clusters(j,3);
        end     
    end
end

RGB_R = reshape(A,size_x,size_y);
RGB_G = reshape(B,size_x,size_y);
RGB_B = reshape(C,size_x,size_y);

image_new = pic;  %新图像
image_new(:,:,1) = RGB_R;
image_new(:,:,2) = RGB_G;
image_new(:,:,3) = RGB_B;
imwrite(image_new,'myimage.jpg','jpg');
imshow(uint8(image_new)),title('result when K = 4');

