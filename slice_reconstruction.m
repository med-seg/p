function slice = slice_reconstruction(sz, centres, score)
% This function determines the probability of each pixel in the original image data (3D array), 
% belonging to positive 'pancreas' via predictions made on the (25x25) patch level by 2 level random forest.

% Inputs: 
% sz – size of the original (3D array) image 
% centres – centre positions of (25x25) patches in the original image: 
% 1) 3 columns for x y and z coordinates, 
% 2) number of rows correspond to the number of extracted patches, 
% score – predicted probability for each (25x25) patch.

% Output: 
% 3D array with restored probabilities

slice = zeros(sz);
for i = 1:length(score)
    patchesX = max(1,centres(i,1)-12):min(centres(i,1)+12,sz(1));
    patchesY = max(1,centres(i,2)-12):min(centres(i,2)+12,sz(2));
    patchesZ = centres(i,3);
    slice(patchesX, patchesY,patchesZ) =...
        max(slice(patchesX,patchesY,patchesZ),score(i));
end

end
