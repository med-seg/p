function h = probabilityMapping(data, pbs, slice, limit)

% This function displays a probability distribution for a 2D image 
% using a selected 'slice' number in a 3D image.

% INPUTS: 
% data - 3D array (matrix) containing image data;
% pbs - corresponding 3D array of probabilities (as predicted by cascade
% random forest);
% slice - selected slice number in a 3D array (image data) starting from 1;
% limit - minimum value of probability to be displayed in colour (default value is 0).

% OUTPUT: 
% handle for displaying image in a figure plot.


%{
%%%% begin: cropping of data to remove excess background %%%
dataCropped = zeros(250,360,50);
for i = 1:50
    [thisBlobsBoundingBox,dataSlice] = cropping2(data,i);
    dataCropped(:,:,i) = dataSlice;
end
data = dataCropped;
%%%% end: cropping of data to remove excess background %%%
%}

% default value
if (nargin < 4)
    limit = 0;
end

% This function draws probability colormap for predicted probabilities
if (size(data,3)>1)
    % read slice of image data
    slice_data = data(:,:,slice);  
end
if(size(pbs,3)>1)
    % read predicted probabilities for slice
    slice_pbs = pbs(:,:,slice);    
else
    slice_pbs = pbs;   
end

p_mask = slice_pbs >= limit;

% apply limit mask for visualization of potential solution
slice_pbs(~p_mask) = 0; 

% count limits for all available probabilities
pb_max = max(slice_pbs(:));

levels = 100;
map = floor(slice_pbs.*(levels-1))+1;
color_pbs = ind2rgb(map,jet(levels));

% use HSV for displaying probability colormap for image slice
hsvColors = rgb2hsv(color_pbs);

% image to HSV
rgbImage = ind2rgb(slice_data,gray(255));
hsvImage = rgb2hsv(rgbImage);

% apply colouring
hcomponent =  hsvColors(:,:,1);
tmp = hsvImage(:,:,1);
tmp(p_mask) = hcomponent(p_mask);
hsvImage(:,:,1) = tmp;
scomponent =  hsvColors(:,:,2);
tmp = hsvImage(:,:,2);
tmp(p_mask) = scomponent(p_mask);
hsvImage(:,:,2) = tmp;

% convert HSV back to RGB
IMAGE = hsv2rgb(hsvImage);  
% draw probability coloured image
h = imshow(IMAGE);
% % colourmap for colorbar
colormap(jet(levels));
%display colobar
c = colorbar;

% set limits for colorbar
if pb_max == 0
pb_max = 0.001;
end

set(c, 'YLim', [0 pb_max]);

end

