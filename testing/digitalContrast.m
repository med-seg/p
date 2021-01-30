% This function performs digital contrast enhancement on a 2D image (slice).
% function enhanceImage = digitalContrast(image, gain, cutoff)

% INPUTS:
% image - 2D image (slice) to be processed
% gain (contrastLevel) - controls the actual contrast
% cutOff (thresh) - normalised greyscale value about which contrast is modified.

% OUTPUT:
% enhanceImage - digitally contrast enhanced 2D image (slice)

function enhanceImage = digitalContrast(image, gain, cutoff)

    if isa(image,'uint8')
	newIm = double(image);
    else 
	newIm = image;
    end
    	
    % rescale image to range [0,1]
    newIm = newIm-min(min(newIm));
    newIm = newIm./max(max(newIm));

    % apply sigmoid function
    enhanceImage =  1./(1 + exp(gain*(cutOff-newIm))); 


end
