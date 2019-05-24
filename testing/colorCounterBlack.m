function [N, colors] = colorCounterBlack(imname)

% This function count number N of foreground colors (black color - background) for image with name 'imname' 
% and foreground colors RGB triplets in range [0..255]

    maxv = max(imname(:));
    map = uint8((double(imname) ./ maxv) .* 255);
    %image(mapped_array);
    %colormap(gray(255));

    %[im, map] = imread(imname);                  % read image with colormap 
    
    im = imname;
    
	if ndims(im) == 3                            % checking image format
		columns = reshape(im, [], 3);
		[unique_colors, ~, ~] = unique(columns, 'rows');   % detect colors triplets

        % use for output foreground colors only (black - background)
        
		colors = unique_colors(2:size(unique_colors,1),:);

        % count number of foreground colors
        
		N = size(unique_colors,1) - 1;
	else
% 		unique_colors = uint8(map * 255); % change range [0..1] to [0..255] for detecting colors triplets
% 		
%         % use for output foreground colors only (black - background)
%         
%         colors = unique_colors(2:size(unique_colors,1),:);
%         
%         % count number of foreground colors
%         
% 		N = size(unique_colors,1) - 1;
        
        mjr_index = unique(im);                 % detect colors of image
        %mjr_index(mjr_index == 0) = [];         % remove black color

        colors = mjr_index;
        
        %colors = uint8(map(mjr_index+1,:)); % change range [0..1] to [0..255] for detecting colors triplets
        N = size(colors,1);                     % count number of color classes of image
	end    

end


