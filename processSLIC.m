function processSLIC(Tagged, slic_size, result_folder, isUpdate)

% This function performs SLIC (Simple Linear Iterative Clustering) algorithm in a 3D array (image data)
% and saves SuperpixelCentres, SuperpixelNumber and Labels to a .mat file

    for i = 1:length(Tagged)
        % path to the original file with image data
        filename = [Tagged(i).path Tagged(i).name '.img'];
        % path to the corresponding file to save superpixel data
        result = [result_folder '/' Tagged(i).name '_slic.mat'];
        result = strrep(result,'//','/');

        % if the update is required, then remove existing files    
        if (isUpdate && exist(result,'file'))
           delete(result);
        end

        %if there is no file with the computed results, then run the computation  
        if (~exist(result,'file'))
            data = readImgFile(filename);
            display(sprintf('SLIC processing file %d of %d: ',i,length(Tagged)));
            disp([filename ' log:']);
            [SuperpixelCenters,SuperpixelNumber,Labels] = slic3D(slic_size,data);
            save(result,'SuperpixelCenters');
            save(result,'SuperpixelNumber','-append');
            save(result,'Labels','-append');
        end
    end
end