function [X,Y] = LoadAllPositivePatches(folder,Tagged,N)

% This function loads provided training set of files (N random positive patches).
% This function is similar to LoadNNegativePatches function but only positive patches are loaded.

% Inputs: 
% folder – directory path of where patch files are saved;
% Tagged - list of .tag (annotated ground-truth) files to use for training;
% N - number of patches to load;

% NOTE: If N is not provided, then all positive patches from all files are loaded. 
% Note that there is possible memory problem in this case.

% Outputs: 
% X - patch features; 
% Y – corresponding labels (0 – non pancreas, 1 – pancreas)

    if (nargin < 3)
        N = inf;
    end
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
            if (isinf(N))
                pportion = nP;
            else
                pportion = portion;
            end
            if(nP>0)
                pSelect = randsample(nP,min(pportion,nP));
                X = [X;PatchFeatures(positiveTagsId(pSelect),:)];
                Y = [Y;PatchTags(positiveTagsId(pSelect))];
            end
        end
        if (isinf(N))
            break;
        end
        loaded = length(Y);
        if(loaded == 0)
            error('No training patch data were found.');
        end
    end
end