function [] = digital_contrast_max_flow()

close all;

inputFolderTest = 'path-to-folder-containing-folders-with-deep-learning-model-predictions-for-a-test-image-file';
[Tagged, notTagged] = listPaths(inputFolderTest);
totalData = Tagged;

for tot = 1:length(totalData)
%for tot = 1:1
    
chr = totalData(tot).path;
[AA,~,~,indexNumber] = sscanf(chr,'%*[^0123456789]%d');

nameOfDirIMG = ([inputFolderTest '\max_flow_' num2str(AA(5))]); 
mkdir(nameOfDirIMG);

filename = [totalData(tot).path totalData(tot).name '.img']
data = readImgFile(filename);


% start finding end, start, thresh, gain
%{

dataCopy = data;
dataCopy = imresize(dataCopy, [384 384]);

meanOfWhole = mean2(dataCopy)
originalImage = dataCopy(:,:,25); 
meanOfImage = mean2(originalImage)
diff = meanOfImage - meanOfWhole
meanOfImage12 = mean2(dataCopy(:,:,12))
meanOfImage37 = mean2(dataCopy(:,:,37))

thresh = (0.0036*meanOfImage) + 0.1886;

takeDG = false;
takeSDG = false;
takeB = false;
takeSB = false;
takeLG = false;
takeL = false;

contrastLevel = -1;

if(diff >= 0.9)
    contrastLevel  = -1;
end

if(diff <= 0.5)
    contrastLevel  = -5;
end

if(diff >= 2.29 && diff <= 2.30)
    contrastLevel  = -5;
end

if(diff >= 3.90 && diff <= 4.35)
    contrastLevel  = -5;
end

if(diff >= 3.20 && diff <= 3.30)
    contrastLevel  = -7;
end

if(diff >= -0.54 && diff <= -0.52)
    contrastLevel  = -7;
end

if(diff >= 1.10 && diff <= 1.20)
    contrastLevel  = -10;
end

if(contrastLevel  == -10)
     takeDG = true;
end

if(contrastLevel  == -7)

    condition1 = meanOfWhole >= 38.4 && meanOfWhole <= 38.5 && meanOfImage >= 37.9 && meanOfImage <= 38.0;
    condition2 = meanOfWhole >= 40.2 && meanOfWhole <= 40.3 && meanOfImage >= 43.5 && meanOfImage <= 43.6;

    if(condition1 == true)
        takeDG = true;
    end
    
    if(condition2 == true)
        takeB = true;
    end
    
end

if(contrastLevel  == -5)
    
    takeDG = true;
    
    condition1 = meanOfWhole >= 38.6 && meanOfWhole <= 38.7 && meanOfImage >= 40.9 && meanOfImage <= 41.0;
    condition2 = meanOfWhole >= 41.5 && meanOfWhole <= 41.6 && meanOfImage >= 41.4 && meanOfImage <= 41.5;

    if(condition1 == true || condition2 == true)
        takeLG = true;
    end
    
    condition3 = meanOfWhole >= 58.3 && meanOfWhole <= 58.4 && meanOfImage >= 62.6 && meanOfImage <= 62.7;
    if(condition3 == true)
        takeB = true;
    end
    
end


if(contrastLevel  == -1)
    
    takeDG = true;
    takeLG = true;
    
    condition1 = meanOfWhole >= 38.7 && meanOfWhole <= 38.8 && meanOfImage >= 40.9 && meanOfImage <= 41.0;
    if(condition1 == true)
        takeLG = false;
    end
    
    condition2 = meanOfWhole >= 40.7 && meanOfWhole <= 40.8 && meanOfImage >= 43.3 && meanOfImage <= 43.4;
    condition3 = meanOfWhole >= 38.0 && meanOfWhole <= 38.1 && meanOfImage >= 43.0 && meanOfImage <= 43.1;
    if(condition2 == true || condition3 == true)
        takeLG = false;
        takeSB = true;
    end
    
    condition4 = meanOfWhole >= 40.8 && meanOfWhole <= 40.9 && meanOfImage >= 42.4 && meanOfImage <= 42.5;
    if(condition4 == true)
        takeB = true;
    end
    
    condition5 = meanOfWhole >= 38.9 && meanOfWhole <= 39.0 && meanOfImage >= 40.5 && meanOfImage <= 40.6;
    condition6 = meanOfWhole >= 41.1 && meanOfWhole <= 41.2 && meanOfImage >= 44.2 && meanOfImage <= 44.3;
    if(condition5 == true || condition6 == true)
        takeDG = false;
        takeSDG = true;
    end
    
end

fprintf('Contrast level is: %d\n', contrastLevel);
fprintf('Taking Black? %d\n', takeB);
fprintf('Taking Dark Grey? %d\n', takeDG);
fprintf('Taking Light Grey? %d\n', takeLG);
fprintf('Taking small dark grey? %d\n', takeSDG);
fprintf('Taking small black? %d\n', takeSB);
fprintf('Taking White? %d\n', takeL);

guideStart = 12;
guideEnd = 38;

aX = abs(diff) >= 1.5 && abs(diff) <= 4.4;
aY = meanOfWhole >= 40.2 && meanOfWhole <= 40.9 && meanOfImage >= 42.4 && meanOfImage <= 43.6;
aZ = meanOfWhole >= 58.3 && meanOfWhole <= 58.4 && meanOfImage >= 62.6 && meanOfImage <= 62.7;

if(aX == true && (aY == true || aZ == true))
    guideStart = 9;
    guideEnd = 32;
    
    if(abs(diff) < 2.5)
        guideEnd = 27;
    end
end


bX = abs(diff) >= 0.1 && abs(diff) <= 3.2;
bY = meanOfWhole >= 37.3 && meanOfWhole <= 38.8 && meanOfImage >= 36.4 && meanOfImage <= 41.0;
bZ = meanOfWhole >= 41.1 && meanOfWhole <= 41.6 && meanOfImage >= 41.4 && meanOfImage <= 44.3;

if(bX == true && (bY == true || bZ == true))
    guideStart = 12;
    guideEnd = 35;

    bXX = meanOfWhole >= 38.4 && meanOfWhole <= 38.8 && meanOfImage >= 37.9 && meanOfImage <= 41.0;
    bYY = meanOfWhole >= 37.3 && meanOfWhole <= 37.4 && meanOfImage >= 36.4 && meanOfImage <= 36.5;
    
    if(bXX == true)
      guideEnd = 27;
    end
    
    if(bYY == true)
      guideEnd = 32;
    end

end

cX = abs(diff) >= 1.60 && abs(diff) <= 2.20;
cAX = abs(diff) >= 2.29 && abs(diff) <= 4.00;
cY = meanOfWhole >= 36.1 && meanOfWhole <= 39.0 && meanOfImage >= 38.6 && meanOfImage <= 41.8;

if((cX == true || cAX == true) && cY == true)
    guideStart = 15;
    guideEnd = 34;
    
    cXX = meanOfWhole >= 36.1 && meanOfWhole <= 36.2 && meanOfImage >= 38.6 && meanOfImage <= 38.7;
    cYY = meanOfWhole >= 38.9 && meanOfWhole <= 39.0 && meanOfImage >= 40.5 && meanOfImage <= 40.6;
    
    if(cXX == true)
        guideStart = 14;
        guideEnd = 38;
    end
    
    if(cYY == true)
        guideStart = 13;
        guideEnd = 41;
    end
    
end


dX = abs(diff) >= 0.40 && abs(diff) <= 1.50;
dAX = abs(diff)>= 1.70 && abs(diff) <= 2.20;
dBX = abs(diff)>= 2.27 && abs(diff) <= 2.28;
dY = abs(diff) >= 5.00;
dZ = meanOfWhole >= 33.4 && meanOfWhole <= 36.8 && meanOfImage >= 34.6 && meanOfImage <= 37.3;
dZZ = meanOfWhole >= 38.0 && meanOfWhole <= 42.3 && meanOfImage >= 38.0 && meanOfImage <= 43.3;

if((dX == true || dY == true || dAX == true || dBX == true)...
        && (dZ == true || dZZ == true))
   
   guideStart = 23;
   guideEnd = 42;
   
   dXX = meanOfWhole >= 36.7 && meanOfWhole <= 36.8 && meanOfImage >= 37.2 && meanOfImage <= 37.3; %27-46
   dYY = meanOfWhole >= 38.0 && meanOfWhole <= 38.1 && meanOfImage >= 43.0 && meanOfImage <= 43.1; %22-38
   dZZZ = meanOfWhole >= 39.8 && meanOfWhole <= 39.9 && meanOfImage >= 39.1 && meanOfImage <= 39.2; %18-43
   
    if(dXX == true)
        guideStart = 27;
        guideEnd = 46;
    end
    
    if(dYY == true)
        guideStart = 22;
        guideEnd = 38;
    end
    
    if(dZZZ == true)
       guideStart = 18;
       guideEnd = 43;
   end
    
end

guideStart
guideEnd

%} 
% end end, start, thresh, gain 

dataCropped = [];
for i = 1:size(data,3)
    [thisBlobsBoundingBox,dataSlice] = cropping2(data,i);
    dataSlice = digitalContrast(uint8(dataSlice),contrastLevel, thresh);
    dataCropped(:,:,i) = dataSlice;
end
%dataCopy = dataCropped;

heightx = size(dataCropped,1);
widthx =  size(dataCropped,2);

nperA = round(0.30*heightx);
nperB = round(0.71*heightx);
nperC = round(0.18*widthx);
nperD = round(0.81*widthx);


dataCropped = dataCropped((nperA-20):(nperB),(nperC+30):nperD,:);
%dataCropped = dataCropped((nperA-30):(nperB),(nperC+30):nperD,:);


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


for limitProb = [0.05 0.15 0.35 0.55 0.75 0.95]
%for limitProb = [0.05]
    
constructImage = zeros(250,360,50);

%resultP = load([mainPath '\segnet180_default.mat']);
resultP = load([mainPath '\segnet180_default_' num2str(limitProb) '.mat']);
%resultP = load([mainPath '\segnet180.mat']);
%resultP = load([mainPath '\segnet180_notsorted_' num2str(limitProb) '.mat']);
predict = resultP.predict;

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

     if(takeL == false)
     oneSlice(~predictionMask)= Lt;
     else
     oneSlice(~predictionMask)= Bk;
     end    
%}
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
%imshow(newOneSlice,[]);
%title(num2str(i));
     
constructImage(:,:,i) = newOneSlice;
     
end

constructImage = imfill(constructImage,'holes');
save([nameOfDirIMG '\constructMaxFlow_' num2str(limitProb) '.mat'],'-v7.3','constructImage');

end

end

end