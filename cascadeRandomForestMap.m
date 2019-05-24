function [maskVolume, maskPancreas] = cascadeRandomForestMap(data,pbs,limit)

% This function creates a probability mapping (using cascade random forest predictions) 
% and a binary mask for every slice in a given image volume (3D array), 
% and creates a binary volume to contain the main pancreas regions of interest.

% INPUTS: 
% data - 3D array (matrix) containing image data;
% pbs - corresponding 3D array of probabilities (as predicted by cascade
% random forest);
% limit - minimum probability.

% OUTPUTS:
% maskVolume - 3D array containing the predicted major pancreas region of
% interest as a binary;
% maskPancreas - 3D array containing the predicted major pancreas region of
% interest;

% set default value
if (nargin < 3)
    limit = 0;
end

% initialise binaryVolume
maskVolume = zeros(size(data));
maskPancreas = zeros(size(data));

for sliceN = 1:size(data,3) % sliceN represents slice number in data (3D array of image data)

% read slice of image data
slice_data = data(:,:,sliceN); 

% read predicted probabilities for sliceN
slice_pbs = pbs(:,:,sliceN);    

% create probability distribution for mask using limit
p_mask = slice_pbs >= limit;
slice_pbs(~p_mask) = 0; 

% create binary mask for sliceN
mask = zeros(size(slice_data));
mask (slice_pbs >= limit) = 1;
mask = imfill(mask,'holes');

% mask out (set to 0) the region in the image slice 
% that does not contain binary 1 (mask)
slice_data(~mask) = 0;

% add mask to maskVolume
maskVolume(:,:,sliceN) = mask;
maskPancreas(:,:,sliceN) = slice_data;

end

end

