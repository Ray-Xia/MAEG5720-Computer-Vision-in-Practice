function [ im_rsz ] = my_resize(im, param,scale)
%my_resize: self defined resized function
%   input args: im--input image; 
%               param--'NN' or 'Bilinear';
%               scale--resize scale
%   output args: im_rsz--resized image; 
img=double(im);
[m,n,p]=size(im);
im_rsz=zeros(floor(scale*m),floor(scale*n),p);
[m_rsz, n_rsz, p] = size(im_rsz);
if(strcmp('NN',param)==1)
    for i=1:m
        for j=1:n
            for k=1:p
                im_rsz(floor(scale*(i-1)+1):floor(scale*i),floor(scale*(j-1)+1):floor(scale*j),k)=img(i,j,k);
            end
        end
    end
elseif(strcmp('Bilinear',param)==1)
    for i=1:m_rsz
        for j=1:n_rsz
            for k=1:p
                f_i=floor(i/scale);c_i=ceil(i/scale);
                f_j=floor(j/scale);c_j=ceil(j/scale);
                if (f_i==0)
                    f_i=1;
                end
                if (f_j==0)
                    f_j=1;
                end
                a00=img(f_i,f_j,k);a01=img(f_i,c_j,k);
                a10=img(c_i,f_j,k);a11=img(c_i,c_j,k);
                dx=i/scale-floor(i/scale);
                dy=j/scale-floor(j/scale);
                im_rsz(i,j,k)=a00*(1-dx)*(1-dy)+a01*(1-dx)*dy+a10*dx*(1-dy)+a11*dx*dy;
            end
        end
    end
else
end
    
im_rsz=(uint8(im_rsz));               
end

