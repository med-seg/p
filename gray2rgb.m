function [imageRGB] = gray2rgb(image)
% This function gives a grayscale image an extra dimension
% in order to use colour within the image.

% INPUT: image - grayscale intensity image.
% OUTPUT: imageRGB - image modified to include 3 colour channels (RGB).

[m n]=size(Image);
rgb=zeros(m,n,3);
rgb(:,:,1)=Image;
rgb(:,:,2)=rgb(:,:,1);
rgb(:,:,3)=rgb(:,:,1);
imageRGB=rgb/255;
end
