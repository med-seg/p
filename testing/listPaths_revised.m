function [WithTag, NoTag] = listPaths_revised(folder)
% This function organises .img files (image data files) into two structures: 
% one with .tags (annotated ground-truth) and no .tag files.

% Input: 
% directory path of main folder containing all other folders

% Outputs:
%  - WithTag - structure with info on .img files, which have related .tag files;
%  - NoTag - structure with info on .img files with no .tag files.

% Fields of output structures:
% - name: name of the .img file without extension;
% - path: path to .img file without filename.

WithTag = [];
NoTag = [];
contents = dir(folder); % dir returns struct with entrie filenames in field 'name'
% counters, to be used as indices when adding elements to output structures
iWT = 0; 
iNT = 0;
for iFld=1:length(contents) % loop through folder
    if  strcmp(contents(iFld).name,'.')
        foldername = [folder '/'];
    else
        if(strcmp(contents(iFld).name,'..'))
            continue;
        else
            foldername = [folder '/' contents(iFld).name '/'];
        end
    end
if isdir(foldername)
    imgList = dir([foldername '*.img']); % find .img files in each folder
    if ~isempty(imgList) % if list of .img files is not empty (else leave for next folder)
        for iName=1:size(imgList,1) % file loop
            filename = imgList(iName).name(1:end-4); % extract filename without extension
            if exist([foldername filename '.tag'], 'file') % if .tag file with same name as .img file exists
                iWT = iWT+1; % index increment (for WithTag structure)
                WithTag(iWT).name = filename; % write structure fields
                WithTag(iWT).path = foldername;
            else % no .tag file (all the same for NoTag structure)
                iNT = iNT+1;
                NoTag(iNT).name = filename;
                NoTag(iNT).path = foldername;
            end 
        end % end file loop
    end
end
end % end folder loop

end

