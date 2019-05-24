function [X,Y] = LoadNRandomPatches(folder,Tagged,N)

% This function loads N random patches from the provided training set of files.

% Inputs: 
% folder – directory path to where patches files are saved 
% Tagged - list of .tag (annotated ground-truth) files to use for training
% N – number of patches to load

% Outputs: 
% X - patch features; 
% Y – corresponding labels (0 – non pancreas, 1 – pancreas).

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
            positiveTagsId = find(PatchTags>0);
            nP = length(positiveTagsId);
            positvePortion = ceil(nP/allTags * portion);
            if(nP>0)
                pSelect = randsample(nP,min(positvePortion,nP));
                X = [X;PatchFeatures(positiveTagsId(pSelect),:)];
                Y = [Y;PatchTags(positiveTagsId(pSelect))];
            end
            nN = allTags - nP;
            if(nN > 0)
                negativeTagsId = find(PatchTags==0);
                nSelect = randsample(nN,min(max(portion - positvePortion,0),nN));
                X = [X;PatchFeatures(negativeTagsId(nSelect),:)];
                Y = [Y;PatchTags(negativeTagsId(nSelect))];
            end
        end
        loaded = length(Y);
        if(loaded == 0)
            error('No training patch data were found.');
        end
    end
    id = randperm(loaded);
    X = X(id,:); Y = Y(id,:);
end