function [min_pix, max_pix, mean_i, sigma_square] =...
    setupPDF(Tagged,result_file, isUpdate)

% This function estimates the main parameters of the intensity of image data amongst all the training data. 
% These parameters include minimum and maximum values, and mean and standard deviation.
% PDF = probability density function

% Inputs: 
% Tagged – list of files from training data; 
% result_file – name of .mat file to store computed results for further use; 
% isUpdate – if this value is set to false, then an existing file will be loaded;
% if the isUpdate value is true, then new recomputations will be performed and saved.

% Outputs: 
% min_pix, max_pix – minimum and maximum intensity values in the training data; 
% mean_i – mean value of intensity;
% sigma_square – standard deviation of intensity.

    if (nargin < 3)
        isUpdate = true;
    end

    if (exist(result_file,'file') && isUpdate)
        delete(result_file);
    else
        if(exist(result_file,'file'))
            load(result_file);
            return;
        end
    end

    min_pix = inf; max_pix = 0;
    N = 0;
    S = 0; S2 = 0;
    for i = 1:length(Tagged)
	
        display(sprintf('PDF settings processing file %d of %d: ',i,length(Tagged)));
		
        % path to the original file with image data
        filename = [Tagged(i).path Tagged(i).name '.img'];
        data = readImgFile(filename);
        
        % compute sum of pixels
        S = S + sum(data(:));
        % compute sum of squared pixels
        S2 = S2 + sum(data(:).^2);
        % compute total number of processed pixels
        N = N + numel(data);

        % estimate minimum and maximum value
        min_level = min(data(:));
        max_level = max(data(:));
        min_pix = min(min_level,min_pix);
        max_pix = max(max_level,max_pix);
    end
    
	% compute mean based on sum and count
    mean_i = S/N;
    % compute standard deviation: sum(x-m)^2/N = m^2 - 2*m*sum(X)/N +
    % sum(x*x)/N
    sigma_square = mean_i^2 + S2/N - 2*mean_i/N * S;

    % save computed data in a .mat file
    save(result_file,'min_pix','-mat');
    save(result_file,'max_pix','-append','-mat');
    save(result_file,'sigma_square','-append','-mat');
    save(result_file,'mean_i','-append','-mat');
end

