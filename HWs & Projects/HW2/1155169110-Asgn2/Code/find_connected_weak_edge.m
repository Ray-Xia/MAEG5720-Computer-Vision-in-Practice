% Find weak edges that are connected to strong edges and set them to 1
function[mag] = find_connected_weak_edge(mag, r, c)
    for i = -2:1:2
        for j = -2:1:2
            if (r+i > 0) && (c+j > 0) && (r+i < size(mag,1)) && ...
                    (c+j < size(mag,2)) % Make sure we are not out of bounds
                if (mag(r+i,c+j) > 0) && (mag(r+i,c+j) < 1) % To find connected weak edge and set it to 1
                    mag(r+i,c+j) = 1;
                    mag = find_connected_weak_edge(mag, r+i, c+j);
                end
            end
        end
    end
end