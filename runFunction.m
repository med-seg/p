function [] = runFunction(varargin)

% This function executes the TRAINING stage:
% 1) train random forest on a training dataset of 3D arrays (image
% volumes) to predict major pancreas regions for organ localisation prior
% to deep learning training;
% 2) train and save a deep learning model using major pancreas regions as predicted by
% the random forest.

% add path to vlfeat binary
run('vlfeat-0.9.19/toolbox/vl_setup');  

% path to folder that stores new results
result_folder = '../results';
% path to folder containing original training experimental data 
% (e.g. MRI, CT and ground-truth annotations)
input_folder  = '../path-to-folder-with-input-data';

if (nargin < 1)
    % default
    isUpdate = false;
else
    % the first argument indicates whether to clear any existing results
    % or otherwise
    isUpdate = varargin{1};
end

% This function checks the data folder (input_folder) to find the (.img) image data 
% with and without corresponding ground-truth (.tag), and returns two
% different structures containing .tag and no .tag file paths
[Tagged, NoTag] = listPaths_revised(input_folder);   

% if the result_folder does not exit,
% then create this folder
if (~exist(result_folder,'dir'))
    mkdir(result_folder)
end

%%% 1. Identify Major (Main) Pancreas Region Phase
 
%%% A) Convert all training image data to superpixels of region (slic size)
%%% 32 ppixels. %%%

% convert all 3D image data that have ground-truth annotations or tags (Tagged) into superpixels of cluster slic_size, 
% and save the new superpixel based 3D images (arrays) to result_folder.
slic_size = 32; % size of slic cluster
processSLIC(Tagged, slic_size, result_folder, isUpdate);
 
%%% B) Generate patches of 25 x 25 from image data. %%%

% pdfMap function computes probability distributions of image data. 
% The probability estimation is performed via Gaussian, 
% using kernel density estimates for the positive (pancreas) and negative
% (pancreas) pixel values.
yPlus = pdfMap(20, Tagged, result_folder, isUpdate);
 
% allData  = [Tagged NoTag];
allData  = Tagged;

% for patch formation, use image data files with .tag files
for i = 1:length(allData)
    fprintf('patches features for file %d of %d: ',i,length(allData));
    % path to the original file with data
    filename = [allData(i).path allData(i).name '.img'];
    
    % path to the corresponding file to save tagging data
    patchFile = [result_folder '/' allData(i).name '_patches.mat'];
    patchFile = strrep(patchFile,'//','/');
      
    % path to .mat file with SLIC results
    slicFile = [result_folder '/' allData(i).name '_slic.mat'];
    slicFile = strrep(slicFile,'//','/');
    
    % path to .tag file with tags (annotated ground-truth)
    tagfile = [allData(i).path allData(i).name '.tag'];
    
    % isUpdate = false;
    % if an update is required, then remove existing files    
    if (isUpdate && exist(patchFile,'file'))
       delete(patchFile);
    end
    
    % if there is no file with the computed results, then run the computation  
    if (~exist(patchFile,'file'))
        data = readImgFile(filename);   % read 3D image data for feature extraction
        slicOutput = load(slicFile);    % load SLIC (superpixel) results
                     
        if (exist(tagfile,'file'))
            % read .tag (ground-truth annotated) file
            [~, tags] = tagRead(tagfile);         
            % extract features for 25x25 2D patches
            [PatchFeatures, PatchCenters, PatchTags] = ...
                patchFeatureExtractor(data, slicOutput.Labels,yPlus,tags);
        else
            % extract features for 25x25 2D patches
            [PatchFeatures, PatchCenters, PatchTags] = ...
                patchFeatureExtractor(data, slicOutput.Labels,yPlus);
        end                    

        % save extracted patches in result_folder
        save(patchFile,'PatchFeatures');
        save(patchFile,'PatchTags','-append');
        save(patchFile,'PatchCenters','-append');   
    end
end

%%% C) Train and save cascade random forest. %%%

% set random forest stage parameters
maxS = 3; % maximum number of cascade layers to use 
T = 50; % maximum number of trees to use on each layer

% isUpdate = false;
% path to save random forest classification results
cascadeFile = [result_folder '/cascadeRandomForest.mat'];
if (isUpdate || ~exist(cascadeFile,'file'))
    % train random forest 
    cascade = cascadeRandomForest(T, maxS, allData, result_folder);
    %save classification results
    %save(cascadeFile, 'cascade', '-v7.3');
    save(cascadeFile, 'cascade');
end

% create folders to store .png files of slices containing
% major pancreas regions and corresponding ground-truth
imgFolder = [result_folder '/img'];
gtFolder = [result_folder '/gt'];
mkdir(imgFolder);
mkdir(gtFolder);

% tracking for converting image slices to .png files after 
% applying random forest prediction 
% to generate major region of interest
mostRecent = 1; 

for i = 1:length(allData)
    %path to the original file with 3D image data
    filename = [allData(i).path allData(i).name '.img'];
    %read 3D image data into a 3D array (data)
    data = readImgFile(filename); 
    %path to save 25x25 patch based probabilities
    pbsname = [result_folder '/' allData(i).name '_pbs.mat'];
    fprintf('Estimating patches probability for %s, file %d of %d',...
            allData(i).name,i,length(allData));
    if (~exist(pbsname,'file') || isUpdate)
        %path to the corresponding file to save tagging data
        patchFile = [result_folder '/' allData(i).name '_patches.mat'];
        patchFile = strrep(patchFile,'//','/');
        patchData = load(patchFile,'PatchCenters','PatchFeatures');

        patchCent = patchData.PatchCenters;       
        testFeatures = patchData.PatchFeatures;

        if(~exist('cascade','var'))
            load(cascadeFile, 'cascade');
        end
        
        [~,score] = cascadePredict(cascade, testFeatures);
        predictedProbabilities = slice_reconstruction(size(data), patchCent, score);

        save(pbsname,'predictedProbabilities');
    else
        load(pbsname);
    end
    
    % map random forest predictions to image data either stored in a folder (as
    % .pngs) or in 3D arrays (matrices), and save to a file.
    maskVolumeName = [result_folder '/' allData(i).name '_majorBinary.mat'];
    majorPancreasName = [result_folder '/' allData(i).name '_majorPancreas.mat'];
    
    [maskVolume, majorPancreas] = cascadeRandomForestMap(data,pbs,limit);
    
    save(maskVolumeName,'maskVolume');
    save(majorPancreasName,'majorPancreas');
        
    %{ 
    % uncomment to visualise sample random forest prediction results %
    
    %generate random test slice (2D image) from data (3D image array)
    testSlice = round(1 + (size(data,3)-1)*rand(1,1));    
    %testSlice = 25;
    
    % set figure handle properties
    f = setFigureProperties; 
    % generate probability mapping for testSlice
    probabilityMapping(data, predictedProbabilities, testSlice, 0.5);
    % path and name of figure (imagename)
    imagename = [result_folder '/' allData(i).name '_rf' num2str(testSlice) '.png'];
    % save figure (imagename) in result_folder
    printFigureWithName(f,imagename); 
    %}      
    
    % extract slices from 3D arrays containing image data and ground-truth
    % and save to two separate folders in result_folder
    
    for kk = 1:size(majorPancreas,3)
    im = majorPancreas(:,:,kk);
    im = uint8(im);
    im = gray2rgb(im);
    
    %printName = [imgFolder '/' 'image_' num2str(kk) '.png'];
    printName = [imgFolder '/' 'image_' num2str(mostRecent) '.png'];
    imwrite(im,printName);

    % path to .tag file with tags (annotated ground-truth)
    filenameTag = [allData(i).path allData(i).name '.tag'];
    [~, groundTruthMRI] = tagRead(filenameTag);
    tag = logical(groundTruthMRI(:,:,kk));

    %printNameGT = [gtFolder '/' 'gt_' num2str(kk) '.png'];
    printNameGT = [gtFolder '/' 'gt_' num2str(mostRecent) '.png'];
    imwrite(tag,printNameGT);

    % update tracking
    mostRecent = mostRecent + 1; 

    end
    
end


% 2. Deep Learning Training Phase

%dataDir = 'path/to/dataFolder';
%imDir = fullfile(dataDir,'folder-containing-training-image-slices');
%pxDir = fullfile(dataDir,'folder-containing-training-ground-truth-for-slices');
imDir = imgFolder;
pxDir = gtFolder;

imds = imageDatastore(imDir,'FileExtensions','.png'); % training images are
classNames = ["pancreas","nonpancreas"];
pixelLabelID = [1 0];
pxds = pixelLabelDatastore(pxDir,classNames,pixelLabelID);

% if image and ground-truth slices are not stored in numerical order
% (depending on file name,i.e. 1,10,100 not 01,02...10,11,...100), then
% sort into numerical order
%{

%ImgsX = fullfile(dataDir,'folder-containing-training-image-slices');
%ImgsXGT = pxDir = fullfile(dataDir,'folder-containing-training-ground-truth-for-slices');
ImgsX = imgFolder; ImgsXGT = gtFolder;
Imgs = dir2(ImgsX); ImgsGT = dir2(ImgsXGT);
CopyImgs = Imgs; CopyGT = ImgsGT;
NumImgs = size(Imgs,1);
tracked = 0;
for o = 1:NumImgs
     for s = 1:NumImgs
         imagename = Imgs(s).name;
         numImage = sscanf(imagename,'image_%d');
         if(numImage == o)
             CopyImgs(o) = Imgs(s);
             CopyGT(o) = ImgsGT (s);
             tracked = tracked + 1;
         end
     end  
end  

CopyImgsCell = cell(NumImgs,1);
CopyImgsCellGT = cell(NumImgs,1);
for t = 1:NumImgs 
    folderName = [CopyImgs(t).folder '\' CopyImgs(t).name];
    CopyImgsCell{t,1} = folderName;
    folderNameGT = [CopyGT(t).folder '\' CopyGT(t).name];
    CopyImgsCellGT{t,1} = folderNameGT;
end

imds.Files = CopyImgsCell;
imds = imageDatastore(CopyImgsCell,'FileExtensions','.png');
pxds = pixelLabelDatastore(CopyImgsCellGT,classNames,pixelLabelID);
%}


% a) configure deep learning architecture

%augmenter = imageDataAugmenter('RandXReflection',true,'RandYReflection',true,'RandRotation', [-45 45],...
%    'RandXTranslation',[-20 20],'RandYTranslation',[-20 20], 'RandXShear',[-10 10], 'RandYShear', [-10 10]);

augmenter = imageDataAugmenter('RandXReflection',true,...
    'RandXTranslation',[-40 40],'RandYTranslation',[-40 40]);

pximds2 = pixelLabelImageSource(imds,pxds,...
    'DataAugmentation',augmenter);

imWidth = 'define';
imHeight = 'define';
imageSize = [imHeight imWidth 3]; % convert possible 8 or 16-bit to to 24-bit
numClasses = numel(classNames); % number of classification classes
lgraph = segnetLayers(imageSize,numClasses,'vgg16');

%figure;
%plot(lgraph);

tbl = countEachLabel(pxds);
frequency = tbl.PixelCount/sum(tbl.PixelCount);
bar(1:numel(classNames),frequency)
xticks(1:numel(classNames)) 
xticklabels(tbl.Name)
xtickangle(45)
ylabel('Frequency')

imageFreq = tbl.PixelCount ./ tbl.ImagePixelCount;
classWeights = median(imageFreq) ./ imageFreq;
%pxLayer = pixelClassificationLayerHS('Name','labels','ClassNames',tbl.Name,'ClassWeights',classWeights)
pxLayer = pixelClassificationLayer('Name','labels','ClassNames',tbl.Name,'ClassWeights',classWeights)
%pxLayer = exampleClassificationSSELayer('Name')

lgraph = removeLayers(lgraph,'pixelLabels');
lgraph = addLayers(lgraph, pxLayer);
%lgraph = addLayers(lgraph, exampleClassificationSSELayer)
lgraph = connectLayers(lgraph,'softmax','labels');

%figure;
%plot(lgraph);


% b) configure deep learning training options

maxEpoch = 'define'; %240; %260; %280 %300% %320 %420;
batchSize = 'define'; %10; %1;

% 'LearnRateSchedule','piecewise',...
options = trainingOptions('sgdm',...
    'Momentum',0.9, ...
    'InitialLearnRate',1e-3, ...
    'LearnRateDropFactor',0.5, ...
    'LearnRateDropPeriod',50, ... 
    'MaxEpochs',maxEpoch, ...
    'MiniBatchSize',batchSize, ...
    'CheckpointPath',tempdir, ...
    'Shuffle','every-epoch', ...
    'VerboseFrequency',2,...
    'Plots','training-progress');
%}

% c) train the deep learning model
[net, info] = trainNetwork(pximds2,lgraph,options);

netFile = [result_folder '/netFile.mat'];
%netFile = strrep(netFile,'//','/');
save(netFile,'net');
save(netFile,'info','-append');

    
end

