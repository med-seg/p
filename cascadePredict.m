function [label,score] = cascadePredict(cascade, X)

% This function predicts the probability of a 25 x 25 patch as positive 'pancreas'.

% Inputs: 
% cascade – model previously 'learned' in the cascadeRandomForest function;
% X – patch features to use for probability prediction.

% Outputs:
% label – assign as positive pancreas for patches with score > 0.5;
% score – probability of the patch belonging to positive pancreas class.


% variable for probability
totalScore = zeros(size(X,1),1); 
% start loop through stages
for i = 1:length(cascade.alpha)
    % read stage forest
    forest = cascade.Forest{i};
    % start looping through trees of forest
    for j = 1:length(forest)
        % calculate sum 
        [~,score] = predict(forest{j},X);
        totalScore = totalScore + cascade.alpha(i)/length(forest).*score(:,strcmp(forest{j}.ClassNames,cascade.positiveName));
    end
end
score = totalScore/sum(cascade.alpha);
label = (score > 0.5);
end


