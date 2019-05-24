function cascade = cascadeRandomForest(T, maxS, Tagged, result_folder)

% This function trains a cascaded random forest model.
% This function calculates and saves predicted probabilities for random
% forest cascade.

% Inputs: 
% T - maximum number of trees to use on each layer; 
% maxS - maximum number of cascade layers to use; 
% Tagged - list of files to use for training; 
% result_folder – directory where patches information is stored.

% Outputs: 
% cascade - model structure consisting of the following components: 
% 	1) Forest  -  cell array of tree arrays learnt at each layer; 
% 	2) alpha - weight coefficients to use for each layer; 
% 	3) positiveClassName - label for the class to compute probability of patches.


% Target TP and TN rates for a stage
 tpLimit = 0.99;
 tnLimit = 0.99;

 % maximum size to load at once
 loadSize = 100000;

%%% Calculate numeric labels for classes %%%

% load initial training set for positive
[positiveX,positiveY] = LoadAllPositivePatches(result_folder,Tagged,loadSize);
[negativeX,negativeY] = LoadNNegativePatches(result_folder,Tagged,loadSize);
stageX = [positiveX;negativeX];
stageY = [positiveY;negativeY];

% load initial test set
[positiveX,positiveY] = LoadAllPositivePatches(result_folder,Tagged,loadSize);
[negativeX,negativeY] = LoadNNegativePatches(result_folder,Tagged,loadSize);
testX = [positiveX;negativeX];
testY = [positiveY;negativeY];

Forest = {};
s = 1;
globalFN = 1;
globalFP = 1;
while (s < maxS && (globalFN > 0.001 || globalFP > 0.1) && globalFN > 0)
    fprintf(['Training random forrest,',... 
        'stage %d, global fnr: %.2f, global fpr: %.2f,'],...
        s,globalFN,globalFP);
    id = randperm(length(stageY));
    stageX = stageX(id,:);
    stageY = stageY(id);
    totalScore = zeros(size(testX,1),1);
    totalScoreT = zeros(size(stageX,1),1);
    ForestC = {};
    for t=1:T
        tree = TreeBagger(T,stageX,stageY,'Method','classification');
        
        % append trained classifier
        ForestC = [ForestC {tree}];
        
        % make estimations of performance for test and training sets
        [~,score] = predict(tree,testX);
        totalScore = totalScore + score(:,strcmp(tree.ClassNames,'1'));
        [~,score] = predict(tree,stageX);
        totalScoreT = totalScoreT + score(:,strcmp(tree.ClassNames,'1'));
        newLT = totalScoreT/t > 0.5;
        TNsetT = (stageY == 0) & (newLT == 0);  % True negative set
        FPsetT = (stageY == 0) & (newLT == 1); % False positive
        FNsetT = (stageY == 1) & (newLT == 0); % False negative
        TPsetT = (stageY == 1) & (newLT == 1);
        TNT = nnz(TNsetT);
        FPT = nnz(FPsetT);
        FNT = nnz(FNsetT);
        TPT = nnz(TPsetT);
        % if goal is achieved for the stage, then stop adding trees
        tprT = TPT/max(TPT+FNT,1);
        tnrT = TNT/max(FPT+TNT,1);
        if (tprT>tpLimit && tnrT>tnLimit)
            break;
        end
    end
    % append models from the stage
    Forest = [Forest {ForestC}];
    
    % remove true negative elements from training set
    stageX(TNsetT,:) = [];
    stageY(TNsetT) = [];
    
    % estimation for this stage on test set
    newL = totalScore/t > 0.5;
    TP = nnz((testY == 1) & (newL == 1)); % True positive
    FP = nnz((testY == 0) & (newL == 1)); % False positive
    FN = nnz((testY == 1) & (newL == 0)); % False negative
    precision = TP/(TP+FP);
    recall = TP/(TP+FN);
    
    % weight for the stage linearly normalized
    alpha(s) = exp(2*precision*recall/(precision + recall))/exp(1);
    
    % estimate global error for overall cascade classifier
    cascade = struct('Forest',{Forest},'alpha',alpha,'positiveName','1');
    l = cascadePredict(cascade,testX);
    FPset = (testY == 0) & (l == 1); % False positive
    FNset = (testY == 1) & (l == 0); % False negative
    globalFN = nnz(FNset)/length(l);
    globalFP = nnz(FPset)/length(l);
    
    % add incorrect values from test set to training set
    stageX = [stageX; testX(FPset,:);testX(FNset,:)];
    stageY = [stageY; testY(FPset);testY(FNset)];
    
    % load new test set
    [positiveX,positiveY] = LoadAllPositivePatches(result_folder,Tagged,loadSize);
    [negativeX,negativeY] = LoadNNegativePatches(result_folder,Tagged,loadSize);
    testX = [positiveX;negativeX];
    testY = [positiveY;negativeY]; 
    
    % go to the next stage
    s = s +1;
end
end

