function [DC] = diceScoreX(tags, pbs)
% This function calculates the dice coefficient for predictive segmentation and ground truth.
% This function returns the diceScore (Sorensen-Dice coefficient) for a 2D (slice) or 3D (whole) image. 

% Inputs:
% tags – 3D / 2D array containing (ground-truth) tags for the 3D image of interest;
% pbs – estimated probabilities for the 3D / 2D image of interest.

% Outputs: 
% DC - Dice coefficient score

A = pbs > 0;

% positive tags
B = tags > 0;           

% intersection of A and B
ACB = A .* B ;  
 
% Dice Coefficient (DC)
denominator = (nnz(A) + nnz(B));
if (denominator > 0)
    DC = (2*nnz(ACB))/(nnz(A) + nnz(B));
else
    DC = 0;
end

end