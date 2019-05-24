close all;

% read .img file containing image data
data = readImageFile('path-to-.img'); 
% read .tag file containing annotations (tags) of organ
[~, tags] = tagRead('path-to-.tag'); 

% map tags (annotations) in red onto a selected slice from image (data)
sliceNumber = 25;
imageHandle = tagDraw(data,tags,sliceNumber); % function that maps (draws) tags onto image slice
title (['Ground-truth annotation for slice ' num2str(sliceNumber)]);
