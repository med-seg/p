function imageHandle = tagDraw(dataImg,dataTags,sliceNumber)
% This function draws ground truth with tags for sliceNumber using a red colour
% This function receives dataImg and corresponding dataTags,
% and displays the ground truth for the selected slice (sliceNumber) of the dataImg

% INPUTS: 
% dataImg: 3D image data (.img); 
% dataTags: corresponding 3D array of tags
% sliceNumber: number of selected slice starting from 1 (to, e.g. 50)

% OUTPUT: handle for displayed image on the figure

% ref: http://www.tech-faq.com/hsv.html
% HSV (Hue, Saturation and Value) – defines a type of colour space.
% hue, saturation and value. ‘Value’ is sometimes substituted with ‘brightness’ 
% and then it is known as HSB.
% In HSV, hue represents colour.
% Saturation indicates the range of grey in the color space. It ranges from 0 to 100%. 
% Sometimes the value is calculated from 0 to 1. 
% When the value is ‘0,’ the color is grey and 
% when the value is ‘1,’ the color is a primary colour (RGB).
% Value is the brightness of the colour and varies with colour saturation. 
% It ranges from 0 to 100%. When the value is ‘0’ the color space will be totally black. 
% With the increase in the value, the colour space brightness up and shows various colours.



%% Extract data for display  %%
sliceData = dataImg(:,:,sliceNumber); % read sliceNumber image
sliceTag = dataTags(:,:,sliceNumber); % read sliceNumber tags
tagMask = logical(sliceTag); % read tag mask

%% Construct indexes from data to RGB then convert RGB to HSV %%

rgbImage = ind2rgb(sliceData,gray(255));
% Convert indexed image to RGB image
% RGB = ind2rgb(X,map) converts the indexed image, X, and the corresponding colormap, map, 
% to the truecolor image, RGB. The indexed image, X, is an m-by-n array of integers. 
% The colormap, map, is a three-column array of values in the range [0,1]. 
% Each row of map is a three-element RGB triplet that specifies the red, green, 
% and blue components of a single color of the colormap. 

hsvImage = rgb2hsv(rgbImage);
% Convert RGB colormap to HSV colormap
% hsv_image = rgb2hsv(rgb_image) converts the RGB image to the equivalent HSV image. 
% rgb_image is an m-by-n-by-3 image array whose three planes contain the red, green, 
% and blue components for the image. HSV is returned as an m-by-n-by-3 
% image array whose three planes contain the hue (H), saturation (S), and value (V)
% components for the image.

% Get HSV representation for red colour (to use for ground truth labelling)
red_hue = hsv2rgb([1.0 0 0]); % double RGB red to HSV
% M = hsv2rgb(H) converts a hue-saturation-value (HSV) colormap to a red-green-blue (RGB) colormap. 
% H is an m-by-3 matrix, where m is the number of colors in the colormap. 
% The columns of H represent hue, saturation, and value, respectively. 
% M is an m-by-3 matrix. Its columns are intensities of red, green, and blue, respectively.

%% Apply red tag colouring %%

% Update hue (H) for the ground truth (red) labelling:
% red hue for masked tag region %
tmp = hsvImage(:,:,1); % extract the hue hsv image part                                      
tmp(tagMask) = red_hue(1); % update tmp (all '1' tagMask elements need to be red)
hsvImage(:,:,1) = tmp; %update hue hsv image part

% Update saturation for masked tag region
tmp = hsvImage(:,:,2); % extract the saturation hsv image part      
tmp(tagMask) = 1; % update tmp (1 meaning primary colour and not grey)
hsvImage(:,:,2) = tmp; % update saturation hsv image part    

% convert HSV to RGB 
convertedImage = hsv2rgb(hsvImage);                                 

% draw image with tag ground truth (red colour)
% The [] lets it scale the min value to 0 and the max value to 255.
imageHandle = imshow(convertedImage,[]);

end