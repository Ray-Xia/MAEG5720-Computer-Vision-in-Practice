% Codes written by QIU Yufu
clc;clear;
pic = imread('peppers.jpg');   %����ͼ��

[size_x,size_y,dimension] = size(pic);   %�õ�ͼ���С

%�������
K = 12;  %����صĸ���

%-------------problem 1--------------%
% clusters = zeros(K, dimension); %�������
% clusters(:,1:3) = randi([1,255],K,3);
%-------------problem 1--------------%

%-------------problem 2--------------%
clusters = zeros(K, dimension+2); %�������
clusters(:,1:3) = randi([1,255],K,3);
clusters(:,4) = randperm(size_x,K);
clusters(:,5) = randperm(size_y,K);
%-------------problem 2--------------%

% ��ͼ��RGB3ͨ���ֽ�
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

sample_num = size(A,1);   %���ص����
dist = zeros(1,K);

c = zeros(sample_num,1);
dist_min = 100000 ;
old_clusters = 1000;  %���ó�ֵ����һ�ε����Ĵ�
old_c = 1000 * ones(sample_num,1);     %���ó�ֵ����һ�ε���������

%-------------problem 2--------------%
% kx = 0;   %x��y�ڷ����е�Ȩ�أ�̫��ᵼ�·ֲ���
% ky = 0;
%-------------problem 2--------------%

iter = 30; % �ٶ�������30��
t = 0;
while(1)
%-------------c-------------%
%     t = t+1;  %(c)
%     if (t == iter)  %��ֵ
%        break;
%     end
%---------------------------%

%------------d(i)--------------%   
%     bb = clusters - old_clusters;
%     bbb = norm(bb);
%     t = t + 1;
%     if (bbb < 1)  %��ֵ
%         break;
%     end
%     old_clusters = clusters;
 %-----------------------------%   
 
 %------------d(ii)--------------%   
    cc = c - old_c;    %���б仯
    ccc = norm(cc);
    t = t + 1;
    if (ccc < 0.5)  %��ֵ
        break;
    end
    old_c = c;
 %-----------------------------%   
 
    for i = 1:sample_num  % ����ÿ�����ص�
        dist_min = 100000;
        for j = 1:K   %����ÿ����cj
            % ��������ģ��
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
    
    %���´��� 
    for j = 1:K   %�����д���
        %��ʼ��
        num = 0;
        A_sum = 0;
        B_sum = 0;
        C_sum = 0;
%------------problem 2-------------%
        X_sum = 0;
        Y_sum = 0;
%------------problem 2-------------%
        for i = 1:sample_num  %���������ص�
            if(c(i) == j)  %�ҵ���Ӧ����µ����ص�
                %�����ص��ֵ׼��
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
        %�����ص��ֵ����ֵ������
        clusters(j,1) = A_sum/num;
        clusters(j,2) = B_sum/num;
        clusters(j,3) = C_sum/num;
    %------------problem 2-------------%
        clusters(j,4) = X_sum/num;
        clusters(j,5) = Y_sum/num;
    %------------problem 2-------------%
    end
end

%�������˾͸���ͼ��
for j = 1:K   %�����д���
    for i = 1:sample_num  %���������ص�
        if(c(i) == j)   %��Ӧ����µ����ص�
           A(i) = clusters(j,1);
           B(i) = clusters(j,2);  %��������ͨ��
           C(i) = clusters(j,3);
        end     
    end
end

RGB_R = reshape(A,size_x,size_y);
RGB_G = reshape(B,size_x,size_y);
RGB_B = reshape(C,size_x,size_y);

image_new = pic;  %��ͼ��
image_new(:,:,1) = RGB_R;
image_new(:,:,2) = RGB_G;
image_new(:,:,3) = RGB_B;
imwrite(image_new,'myimage.jpg','jpg');
imshow(uint8(image_new)),title('result when K = 4');

