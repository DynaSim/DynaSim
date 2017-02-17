

function xp_plotimage (xp,scale)
    % xp must be 1D
    
    rgb = imread(xp.data{1});
    
    if nargin > 1
        rgb = imresize(rgb,scale);
    end
    
    imshow(rgb);
    
end


