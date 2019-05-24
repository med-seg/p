function [] = maxFlowUEF()

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
dataMRI = readImgFile(filename);
filenameTag = [totalData(tot).path totalData(tot).name '.tag'];
[~, groundTruthMRI] = tagRead(filenameTag);
groundTruthMRI = imresize(groundTruthMRI, [260 320]);

contrastLevel = -2;  %1 (2); %-10 (1); %-5(3);
dataCopy = dataMRI;
dataCopy = imresize(dataCopy, [260 320]);
originalImage = dataCopy(:,:,40); %25?
meanOfImage = mean2(originalImage)
thresh = (0.0036*meanOfImage) + 0.1886

%%%%
tagCropped2 = [];
dataUpdate2 = [];
for ii = 1:size(dataMRI,3)
    [BlobsBoundingBox,dataSlice] = cropping4(dataMRI,ii);
    tagSlice = groundTruthMRI(:,:,ii);
    tagSlice = imcrop(tagSlice, BlobsBoundingBox);
    tagSlice = imresize(tagSlice, [170 320]);
    tagCropped2(:,:,ii) = tagSlice;
    dataUpdate2(:,:,ii) = dataSlice;
end
groundTruthMRI = tagCropped2;
dataMRI = dataUpdate2;
%%%%

size(dataMRI)
size(groundTruthMRI)

[boxes labels] = imBoundingBox(logical(groundTruthMRI));
zMax = floor(boxes(6))
zMin = round(boxes(5))

ims = zeros(size(dataMRI));
tags = zeros(size(groundTruthMRI));

for kk = 1:size(dataMRI,3)
im = dataMRI(:,:,kk);
im = digitalContrast(uint8(im),contrastLevel, thresh); %-1 or -3 or -5 or -10 dependng on mean of image slice?
%im = uint8(im);
ims(:,:,kk) = im;
%im = im - min(im(:)); % shift data such that the smallest element of A is 0
%im = im /max(im(:)); % normalize the shifted data to 1 
%im = gray2rgb(im);
%printName = [nameOfDirIMG '/' 'image_' num2str(kk) '.png'];
%imwrite(im,printName);

tag = logical(groundTruthMRI(:,:,kk));
tags(:,:,kk) = tag;
%nameOfDirGT = [totalData(tot).path num2str(numImage) '_gt']
%mkdir(nameOfDirGT);
%printNameGT = ['/home/w1595030/UTest/groundTruthUEF/' 'gt_' num2str(mostRecent) '.png'];
%printNameGT = ['/home/w1595030/testing7/mri_data/dTest/DUN_tests/groundTruthTestCropped/' 'gt_' num2str(mostRecent) '.png'];
%printNameGT = [nameOfDirGT '/' 'gt_' num2str(kk) '.png'];
%imwrite(tag,printNameGT);
%}

end

nii = ims;
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

nlab = 5; %4

ulab(1) = 0;
ulab(2) = 0.16; %0.16
ulab(3) = 0.25; %0.2
ulab(4) = 0.4; %0.4
ulab(5) = 0.53; %0.53

[rows, cols, heights] = size(ur);
vol = rows*cols*heights*nlab;

alpha = ones(rows,cols,heights);
cc = 0.25;    %0.25, 0.005, 0.0025
beta = 1e-4;   %-4
steps = 0.03;  %0.11, 0.03
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


for tt = zMin:zMax

    f1 = figure;

    rotating =  uu(:,:,tt);
    tag = tags(:,:,tt);

[B,L] = bwboundaries(tag);
imshow(rotating,[]);
hold on;
for k = 1:length(B)
   boundary = B{k};
   plot(boundary(:,2), boundary(:,1), 'r', 'LineWidth', 1)
end

    %rotating = gray2rgb(rotating);
    printName = [nameOfDirIMG '/' 'image_' num2str(tt) '.png'];
    saveas(f1,printName);
    %imwrite(rotating, printName);
end


end %end of totData
%}

end %end of function