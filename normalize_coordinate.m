function [nx, ny] = normalize_coordinate(x_point, y_point, xlims, ylims)

% This function returns a pair of the normalized coordinates (x, y)
% based on the width and height of a 2D (image) matrix and 
% the relative position of the pair of coordinates x and y. 
% The size parameters xlims, ylims are provided 
% in form of the vector of length 2 showing initial position and last position.

% INPUTS:
% x-point - relative position of x co-ordinate
% y-point - relative position of y co-ordinate
% xlims - vector contains initial and last position
% ylims - vector contains initial and last position

% OUTPUTS: 
% nx and ny - pair of the normalized coordinates (x, y)

    nx = ((x_point-xlims(1))./(xlims(2) - xlims(1)));
    ny = ((y_point-ylims(1))./(ylims(2) - ylims(1)));
end
