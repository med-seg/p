function [X,Y] = LoadNNegativePatches(folder,Tagged,N)

% This function loads N negative patches from the provided training set of files.
% This function is similar to LoadNRandomPatches function but only negative patches are loaded.

% Inputs: 
% folder – directory path to where patches files are saved 
% Tagged - list of .tag (annotated ground-truth) files to use for training
% N – number of patches to load

% Output: 
% X - patches features; 
% Y – corresponding labels (0 – non pancreas, 1 – pancreas)

    numFiles = length(Tagged);
    X = []; Y = [];
    loaded = 0;
    while(loaded < N)
        rest = N - loaded;
        portion = ceil(rest/numFiles);
        for i = 1:length(Tagged)
            patchFile = [folder '/' Tagged(i).name '_patches.mat'];
            patchFile = strrep(patchFile,'//','/');
            load(patchFile,'PatchFeatures','PatchTags');
            if(PatchTags(1)<0)
                continue;
            end
            allTags = numel(PatchTags);
            negativeTagsId = find(PatchTags==0);
            nN = length(negativeTagsId);
            if(nN > 0)
                negativeTagsId = find(PatchTags==0);
                nSelect = randsample(nN,min(portion,nN));
                X = [X;PatchFeatures(negativeTagsId(nSelect),:)];
                Y = [Y;PatchTags(negativeTagsId(nSelect))];
            end
        end
        loaded = length(Y);
        if(loaded == 0)
            error('No training patch data were found.');
        end
    end
end