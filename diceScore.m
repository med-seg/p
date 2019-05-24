function [DC] = diceScore(tags, pbs, limit)
% This function calculates the dice coefficient for predictive segmentation and ground truth.
% This function returns the diceScore (Sorensen-Dice coefficient) for a 2D (slice) or 3D (whole) image. 

% Inputs:
% tags – 3D /2D array containing (ground-truth) tags for the 3D / 2D image of interest;
% pbs – estimated probabilities for the 3D image of interests
% limit - a minimum probability that corresponds to positive (annotated) tag.

% Output: 
% DC- Dice coefficient score

% apply limit for prediction: only probability higher than or equal to limit is considered
A = pbs >= limit;

% positive (annotated) tags
B = tags > 0;           

ACB = A .* B ;  % intersection of A and B
 
% Dice Coefficient
denominator = (nnz(A) + nnz(B));
if (denominator > 0)
    DC = (2*nnz(ACB))/(nnz(A) + nnz(B)) ;
else
    DC = 0;
end

end

