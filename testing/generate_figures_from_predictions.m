
inputFolderTest = 'path-to-folder-containing-folders-with-deep-learning-model-predictions-for-a-test-image-file';
filePattern = fullfile(inputFolderTest);
totalData = dir2(filePattern); 
listName = {totalData.name}';
newList = [];
checkNumberList = [];
up = 1;
N = length(listName);
for i = 1:N
    fileToCheck = [inputFolderTest '\' char(listName(i,:)) '\dl_pancreas.mat'];
    doesItExist = isfile(fileToCheck);
    if(doesItExist == 1)
        newList{up} = fileToCheck;
        checkNumberList{up} = char(listName(i,:));
        up = up + 1;
    end
end
newList = newList';
checkNumberList = checkNumberList';

for totLIST = 1:length(checkNumberList)

close all;

imageNumber = char(checkNumberList{totLIST});
dataDir = [inputFolderTest '\' imageNumber];
imDir = fullfile(dataDir,'img');
imds = imageDatastore(imDir,'FileExtensions','.png'); 

% if image and ground-truth slices are not stored in numerical order
% (depending on file name,i.e. 1,10,100 not 01,02...10,11,...100), then
% sort into numerical order
%{

ImgsX = fullfile(dataDir,'img');
Imgs = dir2(ImgsX);
CopyImgs = Imgs;
NumImgs = size(Imgs,1);
tracked = 0;
for o = 1:NumImgs
     for s = 1:NumImgs
         imagename = Imgs(s).name;
         numImage = sscanf(imagename,'image_%d');
         if(numImage == o)
             CopyImgs(o) = Imgs(s);
             tracked = tracked + 1;
         end
     end  
end  

CopyImgsCell = cell(NumImgs,1);
for t = 1:NumImgs 
    folderName = [CopyImgs(t).folder '\' CopyImgs(t).name];
    CopyImgsCell{t,1} = folderName;
end

imds.Files = CopyImgsCell;
imds2 = imageDatastore(CopyImgsCell,'FileExtensions','.png');
%}
    
%dataDir = [inputFolderTest '\' imageNumber];
result = load([inputFolderTest '\' imageNumber '\dl_pancreas.mat']);
goldStandard = result.goldStandard;
predict = result.predict;
DC = result.DC;
        
locationNew = ([inputFolderTest '\' imageNumber '\FIG']); 
mkdir(locationNew);
        
p1 = isosurface(goldStandard,0.1);
p2 = patch(p1);
set(p2,'FaceColor','red','EdgeColor','none','FaceAlpha',0.3);
daspect([1,1,1])
view(3); 
camlight; lighting gouraud
p3 = isosurface(predict,0.1);
p4 = patch(p3);
set(p4,'FaceColor','green','EdgeColor','none','FaceAlpha',0.5);
title([imageNumber ' DSC(%): ' num2str(DC)]);

testFile = [locationNew '\' imageNumber];
saveas(gcf,testFile,'fig');
saveas(gcf,testFile,'png');
        
[boxes, ~] = imBoundingBox(goldStandard);
zmin = floor(boxes(5));
zmax = round(boxes(6));
r = round((zmax-zmin).*rand(3,1) + zmin);


for q =[r(1) r(2) r(3)]
    
dataSlice40 = readimage(imds2,(q));    
taggedCropped40 = goldStandard(:,:,q);
newConstruct40 = predict(:,:,q);

% calculate dice score of final segmentation
[DC40] = diceScoreX(taggedCropped40,newConstruct40)

figure;%('units','normalized','outerposition',[0 0 1 1]);
imshow(dataSlice40,[]); 
title([imageNumber ' Slice ' num2str(q)]);
set(gcf,'units','normalized','outerposition',[0 0 1 1]);

testFile = [locationNew '\' imageNumber '_Ori_' num2str(q)];
saveas(gcf,testFile,'png');

figure;%('units','normalized','outerposition',[0 0 1 1]);
imshow(dataSlice40,[]); 
hold on;

boundaries = bwboundaries(taggedCropped40);
numberOfBoundaries = size(boundaries, 1);
for k = 1 : numberOfBoundaries
	thisBoundary = boundaries{k};
	plot(thisBoundary(:,2), thisBoundary(:,1), 'r', 'LineWidth', 1);
end

hold on;

boundariesR = bwboundaries(newConstruct40);
numberOfBoundariesR = size(boundariesR, 1);
for k = 1 : numberOfBoundariesR
	thisBoundaryR = boundariesR{k};
	plot(thisBoundaryR(:,2), thisBoundaryR(:,1), 'g', 'LineWidth', 1);
end

hold off;
title([imageNumber ' ' num2str(q) ' DSC: ' num2str(DC40)]);
set(gcf,'units','normalized','outerposition',[0 0 1 1]);

testFile = [locationNew '\' imageNumber '_slice_' num2str(q)];
saveas(gcf,testFile,'png');

end

end
%%% end of no-post-processing

