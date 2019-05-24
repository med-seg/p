close all;

addpath('/home/w1595030/testing7/mri_code/CMF_ML_v1.0/CMF_ML_v1.0');

cropBoundingBox = [97   165   336   223];

%imageNumber = 13;

FolderToOpen = '/home/w1595030/transferWinSCP3/testing_DCM/';

dirinfo = dir(FolderToOpen);
tf = ismember( {dirinfo.name}, {'.', '..'});
dirinfo(tf) = [];  %remove current and parent directory.
allDataTesting = dirinfo;

%{
FolderToOpen = 'D:\images_CT_Segmentations\testing_DCM\choose20\';
filePattern = fullfile(FolderToOpenIM,'PANCREAS*.*');
allDataTesting = dir(filePattern);
%}

for imageNumber = [10 11 13 15 16 17 19 20 21 22 23 25 27 29]
%for imageNumber = [22]

Ct = [];

fn = [allDataTesting(imageNumber).folder '/' allDataTesting(imageNumber).name];
infos = load(fn);
infos = infos.ims;

ims = [];
ims8 = [];

for kk = 1:size(infos,3)
im = infos(:,:,kk);
im = imcrop(im, cropBoundingBox);
ims(:,:,kk) = im;

im8 = uint8(infos(:,:,kk));
%im8 = imcrop(im8, cropBoundingBox);
ims8(:,:,kk) = im8;

end

middle = floor(size(infos,3)/2);
meanOfWholeCropped = mean2(ims);
originalImageCropped = ims(:,:,middle); 
meanOfImageCropped = mean2(originalImageCropped);
diffCropped = meanOfImageCropped - meanOfWholeCropped;

startGuide = 86;
endGuide = 185;

a = diffCropped > 80.9 && diffCropped <= 82.9;
b = diffCropped > 51.4 && diffCropped <= 53.4;
c = diffCropped > 93.9 && diffCropped <= 95.9;

d = diffCropped > 84.0 && diffCropped <= 86.0;
e = diffCropped > 74.2 && diffCropped <= 79.6;

f = diffCropped > 62.8 && diffCropped <= 63.8;
g = diffCropped > 43.0 && diffCropped <= 44.0;
h = diffCropped > 51.1 && diffCropped <= 51.4;
i = diffCropped > 111.2 && diffCropped <= 113.2;

j = diffCropped > 93.3 && diffCropped <= 93.8;
k = diffCropped > 49.4 && diffCropped <= 51.1;
l = diffCropped > 69.5 && diffCropped <= 72.5;
m = diffCropped > 64.5 && diffCropped <= 66.5;

n = diffCropped > 57.2 && diffCropped <= 59.2;
o = diffCropped > 101.5 && diffCropped <= 103.5;

p = diffCropped > 121.0 && diffCropped <= 123.0;
q = diffCropped > 44.1 && diffCropped <= 44.5;
r = diffCropped > 358.0 && diffCropped <= 360.0;
s = diffCropped > 430.0 && diffCropped <= 432.0;

tt = diffCropped > 74.2 && diffCropped <= 79.4;
u = diffCropped > 79.4 && diffCropped <= 80.9;

if(c) endGuide  = 153;  end
if(o || a || l) endGuide  = 177; end
if(h || n || d) endGuide  = 167; end
if(g || k) endGuide  = 170; end
if(q || b) endGuide  = 188; end
if(j) endGuide  = 216; end
if(s) endGuide  = 219; end
if(i || tt) endGuide  = 183; end
if(m) endGuide  = 142; end
if(r) endGuide  = 202; end
if(f) endGuide  = 159; end
if(p) endGuide  = 139; end
if(u) endGuide  = 162; end

if(a || b || c) startGuide = 80; end
if(d || e) startGuide = 83; end
if(f || g || h || i) startGuide = 86; end
if(j) startGuide = 89; end
if(k) startGuide = 92; end
if(l) startGuide = 96; end
if(m) startGuide = 62; end
if(n || o) startGuide = 70; end
if(p) startGuide = 49; end
if(q) startGuide = 105; end
if(r) startGuide = 121; end
if(s) startGuide = 145; end

disp(['imageNumber = ' num2str(imageNumber)])
disp(['startGuide = ' num2str(startGuide)]);
disp(['endGuide = ' num2str(endGuide)]);
   
unique19 = false;
takeLt = false;
takeLG = true;
takeDLG = false;
takeDG = false;
takeBlk = false;

choice1 = (diffCropped > 50 && diffCropped <= 75) ||...
          (diffCropped > 84 && diffCropped <= 85) ||...
          (diffCropped > 90 && diffCropped <= 145); %LG only

choice2 = (diffCropped > 40 && diffCropped <= 50) ||...
          (diffCropped > 93.9 && diffCropped <= 95.9) ||...
          (diffCropped > 350 && diffCropped <= 390); % DLG only
      
choice3 = (diffCropped > 75 && diffCropped <= 80); % DG + DLG

choice4 = (diffCropped > 57 && diffCropped <= 59) ||...
          (diffCropped > 390 && diffCropped <= 450) ||... 
          (diffCropped > 69.5 && diffCropped <= 71.5); %LG + DLG
      
choice5 = (diffCropped > 80 && diffCropped <= 84); % 19 (DLG+LG unique)

if(choice1)
    %do nothing
end

if(choice2)
    takeLG = false;
    takeDLG = true;
end

if(choice3)
    takeLG = false;
    takeDLG = true;
    takeDG = true;
end

if(choice4)
    takeLG = true;
    takeDLG = true;
    takeDG = false;
end

if(choice5)
    unique19 = true;
    takeLG = true;
    takeDLG = true;
    takeDG = false;
end

disp(['unique19 = ' num2str(unique19)]);
disp(['takeLt = ' num2str(takeLt)]);
disp(['takeBlk = ' num2str(takeBlk)]);
disp(['takeLG = ' num2str(takeLG)]);
disp(['takeDLG = ' num2str(takeDLG)]);
disp(['takeDG = ' num2str(takeDG)]);

aaa = size(infos,1);
bbb = size(infos,2);
ccc = size(infos,3);
%ims = zeros(aaa,bbb,ccc);
%sharpenVolume = zeros(aaa,bbb,ccc);

%{
for pp = 1:size(infos,3)
        
a = im2double(ims(:,:,pp)); %// Read in your image
lap = [-1 -1 -1; -1 8 -1; -1 -1 -1]; %// Change - Centre is now positive
resp = imfilter(a, lap, 'conv'); %// Change

%// Change - Normalize the response image
minR = min(resp(:));
maxR = max(resp(:));
resp = (resp - minR) / (maxR - minR);

%// Change - Adding to original image now
sharpened = a + resp;

%// Change - Normalize the sharpened result
minA = min(sharpened(:));
maxA = max(sharpened(:));
sharpened = (sharpened - minA) / (maxA - minA);

%// Change - Perform linear contrast enhancement
%sharpened = imadjust(sharpened, [60/255 200/255], [0 1]);
sharpened = imadjust(sharpened, [30/255 200/255], [0 1]);
sharpenVolume(:,:,pp) = sharpened;
end
%}

%figure;
%imshow(ims(:,:,153),[]);
    
%figure;
%imshow(infos(:,:,153),[]);

%figure;
%imshow(sharpenVolume(:,:,153),[]);


nii = ims8;
ur = double(nii);
umax = max(max(max(ur)));
umin = min(min(min(ur)));
ur = (ur - umin)/(umax-umin);

%ur = ur(21:120, 70:170, 69:170);

% define the required parameters:
%
%   - nlab: the number of labels or regions. 
%
%   - ulab: ulab(i=1...nlab) gives the nlab labels or image models.
%
%   - alpha: the penalty parameter to the total-variation term.
%       For the case without incorporating image-edge weights, alpha is given
%       by the constant everywhere. For the case with image-edge weights,
%       alpha is given by the pixelwise weight function:
%
%       For example, alpha(x) = b/(1 + a*| nabla f(x)|) where b and a are positive
%       constants and |nabla f(x)| gives the strength of the local gradient.
%
%   - cc: gives the step-size of the augmented Lagrangian method.
%       The optimal range of cc is [0.1, 3].
%
%   - errbound: the error bound for convergence.
%
%   - numIter: the maximum iteration number.
%
%   - steps: the step-size for the graident-projection step to the
%       total-variation function. The optimal range of steps is [0,
%       0.13].
%

nlab = 5; %5

ulab(1) = 0;
ulab(2) = 0.16; %0.16
ulab(3) = 0.2; %0.2
ulab(4) = 0.4; %0.4
ulab(5) = 0.53; %0.53

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
%   - u(x,i=1...nlab) is set to be an initial cut, see below.
%
%   - the source flow field ps(x), see below.
%
%   - the nlab sink flow fields pt(x,i=1...nlab), set to be the specified 
%     legal flows.
%
%   - the spatial flow fiels p(x,i=1...nlab) = (pp1(x,i), pp2(x,i), pp3(x,i)), 
%     set to be zero.

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

for limitX = [0.05 0.15 0.35 0.55 0.75 0.95]
    
aaa = size(infos,1);
bbb = size(infos,2);
ccc = size(infos,3);

constructImage = zeros(aaa,bbb,ccc);

chr = allDataTesting(imageNumber).name;
[AA,~,~,indexNumber] = sscanf(chr,'%*[^0123456789]%d');
AA

resultP = load(['/home/w1595030/transferWinSCP3/CT_segmentation_masks/' num2str(AA) '_aug_240_' num2str(limitX) '.mat']);
predict = resultP.predict;

for i = startGuide:endGuide
    
predictionMask = predict(:,:,i);  
rotating =  uu(:,:,i);

[N, coloursI] = colorCounterBlack(rotating);
[sortedColoursI, ~] = sort(coloursI, 'descend');

Lt = sortedColoursI(1);
    
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

if(N == 5)
   Lt = sortedColoursI(1); %0.53
   LG = sortedColoursI(2); %0.4
   DLG = sortedColoursI(3); %0.2
   DG = sortedColoursI(4);%0.16
   Blk = sortedColoursI(5); %0
end


if(N>=3)
    
%choice 1    
if(unique19 == false && takeLt == false && takeBlk == false...
  && takeLG == true && takeDLG == false && takeDG == false)

    rotating(rotating ~= LG) = Lt;
    
end

% choice 2
if(unique19 == false && takeLt == false && takeBlk == false...
  && takeLG == false && takeDLG == true && takeDG == false)

    rotating(rotating ~= DLG) = Lt;
    
end
 
%choice 3
if(unique19 == false && takeLt == false && takeBlk == false...
  && takeLG == true && takeDLG == true && takeDG == false)

    rotating(rotating == Lt) = Lt;
    rotating(rotating == DG) = Lt;
    rotating(rotating == Blk) = Lt;    
end

%choice 4
if(unique19 == false && takeLt == false && takeBlk == false...
  && takeLG == false && takeDLG == true && takeDG == true)

    rotating(rotating == Lt) = Lt;
    rotating(rotating == LG) = Lt;
    rotating(rotating == Blk) = Lt;    
end

%choice 5
if(unique19 == true)
    
    if(i >= 80 && i<= 114)
    	rotating(rotating == Lt) = Lt;
        rotating(rotating == DG) = Lt;
        rotating(rotating == Blk) = Lt;   
    end
    
    if(i>114)
        rotating(rotating ~= LG) = Lt;
    end
end


end


newTemp = logical(zeros(512,512));
newTemp(165:(165+223),97:(97+336)) = predictionMask;
rotating(~newTemp) = Lt; %10,2,3
%rotating(~mask) = Blk; %1

%figure;
%imshow(rotating,[]);
%title(num2str(i));

rotatingBlobs = logical(rotating);
ccc = bwconncomp(rotatingBlobs, 8);
numBlobs = ccc.NumObjects;

if(numBlobs >=1)
    
newOneSlice = rotating;
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

%%%%% get Background and flip if necessary
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
%%%%% get Background and flip if necessary

%figure;
%imshow(newOneSlice);
%title(['New slice AFTER Binary' num2str(i)]);

constructImage(:,:,i) = newOneSlice;

else
    fprintf('No segmentation in this slice.\n');
    newOneSlice = zeros(aaa,bbb);
    constructImage(:,:,i) = newOneSlice;
end

end

constructImage = imfill(constructImage,'holes');

fileX = ['/home/w1595030/transferWinSCP3/newCT/constructImageOriginal_CT' num2str(imageNumber) '_aug_240_' num2str(limitX) '.mat'];
save(fileX,'-v7.3','constructImage');

end

end
