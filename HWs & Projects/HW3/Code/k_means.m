function [img_change,label_change]=k_means(img,k,iterations)
%% randomly initialize the centre
[Row,Col,r]=size(img);
img=double(img);
random_Row=randperm(Row);
random_Col=randperm(Col);
c=zeros(k,1,r);
label=zeros(Row,Col);
for i=1:k
    c(i,1,:)=img(random_Row(i),random_Col(i),:);
end
iteration=1;
while iteration<iterations
   iteration=iteration+1; 
%% compute the distance
   dist_min=100000;
    for i=1:Row
        for j=1:Col
            for kk=1:k
                dist=sqrt((img(i,j,1)-c(kk,1,1))^2+(img(i,j,2)-c(kk,1,2))^2+(img(i,j,3)-c(kk,1,3))^2);
                if dist<dist_min
                    label(i,j)=kk;
                    dist_min=dist;
                end
            end
            dist_min=100000;
        end
    end
%% compute the new centre
    add_r=0;
    add_g=0;
    add_b=0;
    total_number=0;
    for kk=1:k
        for i=1:Row
            for j=1:Col
                if label(i,j)==kk
                    add_r=add_r+img(i,j,1);
                    add_g=add_g+img(i,j,2);
                    add_b=add_b+img(i,j,3);
                    total_number=total_number+1;
                end
            end
        end
        c(kk,1,1)=add_r/total_number;
        c(kk,1,2)=add_g/total_number;
        c(kk,1,3)=add_b/total_number;
        total_number=0;
        add_r=0;
        add_g=0;
        add_b=0;
    end
end
for kk=1:k
        for i=1:Row
            for j=1:Col
                if label(i,j)==kk
                    img(i,j,:)=c(kk,1,:);
                end
            end
        end
end
img_change=img;
label_change=label;

end