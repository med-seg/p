function [] = digital_contrast_max_flow_limits()

close all;

inputFolderTest = 'path-to-folder-containing-folders-with-deep-learning-model-predictions-for-a-test-image-file';
[Tagged, notTagged] = listPaths_revised(inputFolderTest);
totalData = Tagged;

for tot = 1:length(totalData)
%for tot = 1:1

% if necessary, modify two lines below
% to extract a string or number that denotes the test image volume    
chr = totalData(tot).path;
[AA,~,~,indexNumber] = sscanf(chr,'%*[^0123456789]%d');

% create folder to store segmentation results
nameOfDirIMG = ([inputFolderTest '\max_flow_' num2str(AA(5))]); 
mkdir(nameOfDirIMG);

% read image data (3D array)
filename = [totalData(tot).path totalData(tot).name '.img']
data = readImgFile(filename);

% set default values for 
% digital contrast enhancement

contrastLevel = -1;
thresh = 0.5;

% crop data to remove background (black)
dataCropped = [];
for i = 1:size(data,3)
    [thisBlobsBoundingBox,dataSlice] = cropping2(data,i);
    dataCropped(:,:,i) = dataSlice;
end

% compute mean input values for artificial neural network 
% that predicts the gain (contrastLevel) and threshold (thresh) 
% for digital contrast enhancement (DCE).

m1 = mean2(dataCropped);

middle = round(0.5*size(data,3));
m2 = dataCropped(:,:,middle);

quarterSlice = round(0.25*size(data,3));
m3 = dataCropped(:,:,quarterSlice);

threeQuarterSlice = round(0.75*size(data,3));
m4 = dataCropped(:,:,threeQuarterSlice);

m5 = dataCropped(:,:,1:quarterSlice);
m5 = dataCropped(:,:,threeQuarterSlice:size(data,3));

% run artificial neural network to predict 
% the gain (contrastLevel) and threshold (thresh) OR
% run artificial neural network to predict the gain (contrastLevel) 
% and then run the linear regression to predict threshold (thresh)

result_gtl = load('path-to-where-artificial-neural-network-is-stored.mat');
net = result_gtl.net;
output = net([all_six_inputs]);
contrastLevel = output(1);
thresh = output(2);

% apply digital contrast
for i = 1:size(dataCropped,3)
    %[thisBlobsBoundingBox,dataSlice] = cropping2(data,i);
    dataSlice = digitalContrast(uint8(dataSlice),contrastLevel, thresh);
    dataCropped(:,:,i) = dataSlice;
end


%guideStart = 12;
%guideEnd = 38;

%{

%crop image to remove non-pancreas background
heightx = size(dataCropped,1);
widthx =  size(dataCropped,2);

nperA = round(0.30*heightx);
nperB = round(0.71*heightx);
nperC = round(0.18*widthx);
nperD = round(0.81*widthx);

dataCropped = dataCropped((nperA-20):(nperB),(nperC+30):nperD,:);
%dataCropped = dataCropped((nperA-30):(nperB),(nperC+30):nperD,:);

%}

nii = dataCropped;
ur = double(nii);
umax = max(max(max(ur)));
umin = min(min(min(ur)));
ur = (ur - umin)/(umax-umin);

nlab = 4; %4

ulab(1) = 0;
ulab(2) = 0.16; %0.16
ulab(3) = 0.4; %0.4
ulab(4) = 0.53; %0.53

takeDG = false;
takeSDG = false;
takeB = false;
takeSB = false;
takeLG = false;
takeL = false;

% predict which colours (ulab - labels) to keep (pancreas) and which to disgard (non-pancreas)
% this additional stage willadd to the final accuracy following masking of the VBM (volume binary mask)
% that represents the predicted pancreas mask as a volume.

% run artificial neural network model that makes prediction
% based on mean of image data (3D array) and difference in grayscale
% intensities
% update take[DG, SDG, B, SB, LG, L]

[rows, cols, heights] = size(ur);
vol = rows*cols*heights*nlab;

alpha = ones(rows,cols,heights);
cc = 0.25;    %0.25, 0.005, 0.0025
beta = 1e-4;   %-4
steps = 0.03;  %0.11
iterNum = 400; %200, latest 400

% build up the data terms
for i=1:nlab
    Ct(:,:,:,i) = abs(ur - ulab(i))*10;
end

% set the initial values:
u = zeros(rows,cols,heights,nlab);
pt = zeros(rows,cols,heights,nlab);

[ps,I] = min(Ct, [], 4);

for k=1:rows
    for j=1:cols
        for l=1:heights
            pt(k,j,l,:) = ps(k,j,l);
            u(k,j,l,I(k,j,l)) = 1;
        end
    end
end

divp = zeros(rows,cols,heights,nlab);

pp1 = zeros(rows, cols+1, heights, nlab);
pp2 = zeros(rows+1, cols, heights, nlab);
pp3 = zeros(rows, cols, heights+1, nlab);

erriter = zeros(iterNum,1);

tic
for i = 1:iterNum
    
    pd = zeros(rows,cols,heights);
    
    % update the flow fields within each layer i=1...nlab
    
    for k= 1:nlab
        
        % update the spatial flow field p(x,i) = (pp1(x,i), pp2(x,i), pp3(x,i)):
        % the following steps are the gradient descent step with steps as the
        % step-size.
        
        ud = divp(:,:,:,k) - (ps - pt(:,:,:,k) + u(:,:,:,k)/cc);
        
        % update the component pp1(x,i)
        
        pp1(:,2:cols,:,k) = steps*(ud(:,2:cols,:) - ud(:,1:cols-1,:)) + pp1(:,2:cols,:,k); 
        
        % update the component pp2(x,i)
        
        pp2(2:rows,:,:,k) = steps*(ud(2:rows,:,:) - ud(1:rows-1,:,:)) + pp2(2:rows,:,:,k);
        
        % update the component pp3(x,i)
        
        pp3(:,:,2:heights,k) = steps*(ud(:,:,2:heights) - ud(:,:,1:heights-1)) + pp3(:,:,2:heights,k);
        
        % the following steps give the projection to make |p(x,i)| <= alpha(x)
        
        gk = sqrt((pp1(:,1:cols,:,k).^2 + pp1(:,2:cols+1,:,k).^2 + pp2(1:rows,:,:,k).^2 ...
            + pp2(2:rows+1,:,:,k).^2 + pp3(:,:,1:heights,k).^2 + pp3(:,:,2:heights+1,k).^2)*0.5);

        gk = double(gk <= alpha) + double(~(gk <= alpha)).*(gk ./ alpha);
        gk = 1 ./ gk;
        
        % update the component pp1(x,i)
        
        pp1(:,2:cols,:,k) = (0.5*(gk(:,2:cols,:) + gk(:,1:cols-1,:))).*pp1(:,2:cols,:,k); 
        
        % update the component pp2(x,i)
        
        pp2(2:rows,:,:,k) = (0.5*(gk(2:rows,:,:) + gk(1:rows-1,:,:))).*pp2(2:rows,:,:,k);
        
        % update the component pp3(x,i)
        
        pp3(:,:,2:heights,k) = (0.5*(gk(:,:,2:heights) + gk(:,:,1:heights-1))).*pp3(:,:,2:heights,k);
        
        % recompute the divergence field divp(x,i)
        
        divp(:,:,:,k) = pp1(:,2:cols+1,:,k)-pp1(:,1:cols,:,k)+pp2(2:rows+1,:,:,k)-...
            pp2(1:rows,:,:,k) + pp3(:,:,2:heights+1,k) - pp3(:,:,1:heights,k);
        
        % update the sink flow field pt(x,i)
        
        ud = - divp(:,:,:,k) + ps + u(:,:,:,k)/cc;
        pt(:,:,:,k) = min(ud, Ct(:,:,:,k));
        
        % pd: the sum-up field for the computation of the source flow field
        %      ps(x)
        
        pd = pd + (divp(:,:,:,k) + pt(:,:,:,k) - u(:,:,:,k)/cc);
        
    end
    
    % updata the source flow ps

    ps = pd / nlab + 1 / (cc*nlab);
    
	% update the multiplier u
    
    erru_sum = 0;
    for k = 1:nlab
	    erru = cc*(divp(:,:,:,k) + pt(:,:,:,k) - ps);
	    u(:,:,:,k) = u(:,:,:,k) - erru;
        erru_sum = erru_sum + sum(sum(sum(abs(erru))));
    end
    
    % evaluate the avarage error  
    erriter(i) = erru_sum/vol;
    
    if erriter(i) < beta
        break;
    end
    
end

toc
timet = toc

msg = sprintf('number of iterations = %u. \n', i);
disp(msg);

[um,I] = max(u, [], 4);

% reconstructing the image with the computed label-indicator function I(x)
uu = zeros(rows,cols,heights);

for k=1:rows
    for j=1:cols
        for l=1:heights
            uu(k,j,l) = ulab(I(k,j,l));
        end
    end
end

thresholds = linspace(0.05,0.95,10);
for yy = 1:length(thresholds)

limitProb = thresholds(yy);
%for limitProb = [0.05 0.15 0.35 0.55 0.75 0.95]
%for limitProb = [0.05]
    
constructImage = zeros(size(data));

resultP = load(['path-where-deep-learning-prediction-is-stored' '\dl_pancreas_limit' num2str(limitProb) '.mat']);
predict = resultP.predict;

[boxes, ~] = imBoundingBox(logical(predict));
guideStart = round(boxes(5));
guideEnd = round(boxes(6));

%for i = 1:size(data,3) 
for i = guideStart:guideEnd   
    
     predictionMask = predict(:,:,i);
     oneSlice = uu(:,:,i);
     
     Lt = max(oneSlice(:));
     Bk = min(oneSlice(:));
     
    [N, coloursI] = colorCounterBlack(oneSlice);
    [sortedColoursI, ~] = sort(coloursI, 'descend');
    
     if(N == 2)
        Lt = sortedColoursI(1);
        Blk = sortedColoursI(2);
     end
     
     if(N == 3)
        Lt = sortedColoursI(1);
        LG = sortedColoursI(2);
        DG = sortedColoursI(3);
        Blk = sortedColoursI(3);
     end
     
     if(N == 4)
        Lt = sortedColoursI(1);
        LG = sortedColoursI(2);
        DG = sortedColoursI(3);
        Blk = sortedColoursI(4);
     end
     
    
    if(N==4 && takeL == true)
         oneSlice(oneSlice ~= Lt) = Blk;
     end
     
     if(N==4 && takeL == false && takeLG == false && takeDG == true && takeB == true)
         oneSlice(oneSlice == LG) = Lt;
     end

     if(N==4 && takeL == false && takeLG == false && takeDG == false && takeB == true...
             && takeSDG == false)
         oneSlice(oneSlice ~= Blk) = Lt; 
     end
     
     if(N==4 && takeL == false && takeLG == false && takeSDG == true && takeB == true...
             && takeDG == false)
         oneSlice(oneSlice == LG) = Lt;
     end
     
     if(N==4 && takeL == false && takeLG == false && takeDG == false && takeB == true...
             && takeSDG == true)
         oneSlice(oneSlice == LG) = Lt;
     end
     
     if(N==4 && takeL == false &&...
             takeB == false && takeLG == false && takeDG == true && takeSB == true)
         oneSlice(oneSlice == LG) = Lt;
     end
     
     if(N==4 && takeL == false && takeLG == false && takeDG == true && takeB == false...
             && takeSB == false)
         oneSlice(oneSlice ~= DG) = Lt;
     end
     
     if(N==4 && takeL == false &&...
             takeLG == true && takeDG == true && takeB == false && takeSB == false)
         oneSlice(oneSlice == Blk) = Lt;
     end
     
     if(N==4 && takeL == false &&...
             takeLG == true && takeDG == false && takeB == false &&...
             takeSB == false && takeSDG == true)
         oneSlice(oneSlice == Blk) = Lt;
     end
  
     %%%%
     
     if(N==3 && takeL == true)
         oneSlice(oneSlice ~= Lt) = Blk;
     end

     if(N==3 && takeL == false && takeLG == true && takeDG == true && takeB == false && takeSB == false)
         % do nothing!
     end
     
     if(N==3 && takeL == false &&...
             takeLG == false && takeDG == true && takeB == false && takeSB == false)
         oneSlice(oneSlice ~= DG) = Lt;
     end
     
     if(N==3 && takeL == false &&...
             takeLG == false && takeDG == false && takeB == true && takeSDG == true)
         % do nothing!
     end
     
     if(N==3 && takeL == false && takeLG == false && takeDG == false && takeB == true)
         oneSlice(oneSlice ~= Blk) = Lt; 
     end
     
     if(N==3 && takeL == false && takeLG == false && takeDG == true && takeB == true)
         oneSlice(oneSlice == LG) = Lt;
     end
     
     if(N==3 && takeL == false &&...
             takeB == false && takeLG == false && takeDG == true && takeSB == true)
         oneSlice(oneSlice == LG) = Lt;
     end
     
     %}
    
     %{
     if(takeL == false)
     oneSlice(~predictionMask)= Lt;
     else
     oneSlice(~predictionMask)= Bk;
     end 
     %}
     
	%{
	if(takeL == false)
	   oneSliceBig = ones(250,360);
	   oneSliceBig(:)= Lt; 
	   oneSliceBig((nperA-20):(nperB),(nperC+30):nperD) = oneSlice;
	   oneSlice = oneSliceBig;
	else
		oneSliceBig = zeros(250,360);
		oneSliceBig(:)= Bk;    
		oneSliceBig((nperA-20):(nperB),(nperC+30):nperD) = oneSlice;
		oneSlice =  oneSliceBig;
	end
	%}     

     if(takeL == false)
     oneSlice(~predictionMask)= Lt;
     else
     oneSlice(~predictionMask)= Bk;
     end    
    
	 %figure;
     %imshow(oneSlice,[]);
     %title(num2str(i));
     
     newOneSlice = oneSlice;
     newOneSlice = imcomplement(newOneSlice); %10,2,3
    
%figure;
%imshow(newOneSlice,[]);
%title(['Before slice' num2str(i)]);

darkest =  min(newOneSlice(:));
for ii = 1:size(newOneSlice,1)
    for jj = 1:size(newOneSlice,2)
        pixel = newOneSlice(ii,jj);
        if(pixel ~= darkest) 
            newOneSlice(ii,jj) = 1; 
        end
     end
end


%figure;
%imshow(newOneSlice,[]);
%title(['New slice before Binary' num2str(i)]);

newOneSlice = imbinarize(double(newOneSlice));
newOneSlice = imfill(newOneSlice,'holes');
newOneSlice = bwareaopen(newOneSlice,10); 

%figure;
%imshow(newOneSlice,[]);
%title(['New slice before Binary2' num2str(i)]);

%%%%% start: get background and flip if necessary

leftCol = newOneSlice(:,1); 
numLC = numel(leftCol); %e.g 250
sumLC = sum(leftCol);   %e.g. 250
ratioLC = sumLC/numLC;  %e.g. 1 == white

rightCol = newOneSlice(:,end);  
numRC = numel(rightCol);
sumRC = sum(rightCol);
ratioRC = sumRC/numRC;

topRow = newOneSlice(1,:); 
numTR = numel(topRow);
sumTR = sum(topRow);
ratioTR = sumTR/numTR;

bottomRow = newOneSlice(end,:);
numBR = numel(bottomRow);
sumBR = sum(bottomRow);
ratioBR = sumBR/numBR;

%over 90% are white, then flip
if(ratioLC >= 0.9 || ratioRC >= 0.9 || ratioTR >= 0.9 || ratioBR >= 0.9) 
    fprintf('Flip\n');
    newOneSlice = ~newOneSlice;
else
    %do nothing
end

%%%%% end: get background and flip if necessary

%figure;
%imshow(newOneSlice,[]);
%title(num2str(i));
     
constructImage(:,:,i) = newOneSlice;
     
end

constructImage = imfill(constructImage,'holes');
save([nameOfDirIMG '\constructMaxFlow_' num2str(limitProb) '.mat'],'-v7.3','constructImage');

end

end

end