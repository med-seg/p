function [mod_im, left, right, bottom, up] = table_remove(im)

% This function removes excess 'border' pixels (pixels that belong to non-target object background).
% This function cuts the meaningful part of the 2D image (single slice in 3D image data). 
% The estimation is performed on the average values of the intensity of the image: 
% 	it is assumed that the regions with no intensity values higher than 0.1 of overall mean can be removed. 
% 	Small pixel regions which usually correspond to noise are removed. 
% 	The minimum and maximum filtered pixel coordinate in each dimension are selected and used as the new border of the image. 

% INPUT: 
% im - initial (2D) image

% OUTPUTS: 
% mod_im - cropped (cut) image;
% left - left border of cutting;
% right - right border of cutting;
% bottom - bottom border of cutting;
% up - up border of cutting.

% binary mask for image
bw = double(im>(mean(im(:))*0.1));

% remove tiny objects that usually correspond to noise
bw2 = bwareaopen(bw, 10);

% calculate borders for table removal
x_borders = any(bw2,1);   % borders of X dimension
y_borders = any(bw2,2);   % borders of Y dimension

left = find(any(x_borders,1),1,'first'); % left border
right = find(any(x_borders,1),1,'last'); % right border

bottom = find(any(y_borders,2),1,'first'); % bottom border
up = find(any(y_borders,2),1,'last'); % up border

% apply borders to image and generate cropped image
mod_im = im(bottom:up,left:right);

end

