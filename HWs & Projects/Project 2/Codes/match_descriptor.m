function [matches,scores] = match_descriptor(d1,d2)

% Initialization
thresh = 1.5; 
num_of_d1 = size(d1,2);% number of descriptors in image1
num_of_d2 = size(d2,2);% number of descriptors in image2
fewerDescriptors = min(num_of_d2,num_of_d1);
moreDescriptors = max(num_of_d2,num_of_d1); % maximum number of matches
scores = zeros(1,fewerDescriptors); 
matches = zeros(2,fewerDescriptors);
num_of_matches = 0;
temp = zeros(2,moreDescriptors); 

for i = 1:num_of_d1 
    for j = 1:num_of_d2 
        temp(1 , j) = j;
        temp(2 , j) = norm(double(d1(:,i)-d2(:,j)),2);     
    end

    temp = temp(:,1:fewerDescriptors);    
    sorted = transpose(sortrows(transpose(temp),2));% sort 2nd column of temp transpose
    min1 = sorted(:,1); % 2nd row is the minimum distance, 1st row is descriptor index
    min2 = sorted(:,2); % the second minimum distance
    NNDR = min2(2,:) / min1(2,:);

    if NNDR >= thresh        
        num_of_matches = num_of_matches +1;
        matches(1,num_of_matches) = i;
        matches(2,num_of_matches) = min1(1,:);
        scores(:,num_of_matches) = norm(double(d1(:,i) - d2(:,min1(1,:))),2);

    end
end

matches = matches(:,1:num_of_matches);

end