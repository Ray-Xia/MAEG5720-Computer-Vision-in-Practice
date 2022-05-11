function output_image = box_filter(input_image,filter)

[r,c]=size(input_image);
[m,n]=size(filter);
h=rot90(filter,2);
center=floor((size(h)+1) / 2); %floor向下取整
left=center(2)-1;
right=n-center(2);
top = center(1) - 1;
bottom = m - center(1);
Rep = zeros(r + top + bottom, c + left + right);

for x = 1 + top : r + top
    for y = 1 + left : c + left
        Rep(x,y) = input_image(x - top, y - left);
    end
end

output_image = zeros(r , c);
for x = 1 : r
    for y = 1 : c
        for i = 1 : m
            for j = 1 : n
                q = x - 1;
                w = y -1;
                output_image(x, y) = output_image(x, y) + (Rep(i + q, j + w) * h(i, j));
            end
        end
    end
end

end