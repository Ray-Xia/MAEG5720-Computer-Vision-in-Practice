% The following codes are from
% https://github.com/guyufei96/Harris-Corner-Detector-K-Means-Eigenface/blob/master/my_kmeans_function.m
function [ Ind,  centroids ] =my_kmeans_function(dataSet, k)  
    [m,n] = size(dataSet);      
    centroids = zeros(k, n);  
    for j = 1:n  
        minJ = min(dataSet(:,j));  
        rangeJ = max(dataSet(:,j))-min(dataSet(:,j));  
        centroids(:,j) = minJ+rand(k,1)*rangeJ;
    end  
	  
    subCenter = zeros(m,2);  
    Ind = zeros(m, 1);
    change = 1;
    while change == 1  
        change = 0;  
        for i = 1:m  
            minDist = inf;  
            minIndex = 0;  
            for j = 1:k  
                 x = dataSet(i,:);
                 y = centroids(j,:);
                 dist= (x-y)*(x-y)';  
                 if dist < minDist  
                     minDist = dist;  
                     minIndex = j;  
                 end  
            end  
            if subCenter(i,1) ~= minIndex  
                change = 1;  
                subCenter(i,:)=[minIndex, minDist]; 
                Ind = subCenter(:,1);
            end          
        end  
       
        for j = 1:k  
            sum = zeros(1,n);  
            r = 0; %ÊýÁ¿  
            for i = 1:m  
                if subCenter(i,1) == j  
                    sum = sum + dataSet(i,:);  
                    r = r+1;  
                end  
            end  
            centroids(j,:) = sum./r;  
        end 

    end  

end
