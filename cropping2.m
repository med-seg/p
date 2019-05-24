function [thisBlobsBoundingBox,subImage] = cropping2(data, sliceNumber)

% This function crops a chosen 2D image (slice) for given 3D image data
% and returns the cropped 2D image and cropping bounding box (min-x, min-y, width, height).

% Inputs:
% data - original 3D array (matrix) of image data;
% sliceNumber - selected slice number (z-axis) that represents a 2D image (axial slice) in the 3D image data.

% Outputs:
% subImage - cropped 2D image (slice)
% thisBlobsBoundingBox - measurements needed to crop original 2D image (originalImage) into subImage


data = imresize(data, [384 384]); % resize image width and height to 384 x 384
originalImage = data(:,:,sliceNumber);
thresholdValue = 20;
binaryImage = originalImage > thresholdValue;
binaryImage = imfill(binaryImage, 'holes');
[labeledImage, numberOfBlobs] = bwlabel(binaryImage);
blobMeasurements = regionprops(labeledImage, 'area');
allAreas = [blobMeasurements.Area];
[sortedAreas, sortIndexes] = sort(allAreas, 'descend');
blobMeasurements = regionprops(labeledImage, originalImage, 'all');
thisBlobsBoundingBox = blobMeasurements(sortIndexes(1:1)).BoundingBox;
subImage = imcrop(originalImage, thisBlobsBoundingBox);
subImage = imresize(subImage, [250 360]);

end
