function [C, K, l] = slic3D(regionSize, data)

% This function performs SLIC (Simple Linear Iterative Clustering) algorithm in a 3D image (3D array)

% Inputs:
% regionSize – size of superpixel; 
% data – 3D array containing data loaded from .img (image data) file

% Outputs:
% C = SupepixelCentres: superpixels centres (average color, x, y and z coordinates) for the clusters found with SLIC algorithm;
% K = SupepixelNumber: number of superpixels;
% l = Labels: corresponding labels for each pixel of the input data. Each label corresponds to the 
% particular centre of the cluster found with the SLIC algorithm.

% Outputs:
% C – array of centres for segments with size 4xK; 
% K – number of extracted superpixels; 
% l – labels for each pixel in 3D array data (label is a number that corresponds centre in C)

%{
%%%% begin: cropping of data to remove excess background %%%
dataCropped = zeros(250,360,50);
for i = 1:50
    [thisBlobsBoundingBox,dataSlice] = cropping2(data,i);
    dataCropped(:,:,i) = dataSlice;
end
data = dataCropped;
%%%% end: cropping of data to remove excess background %%%

%}

% determine dimensions of the image data
[width, height, depth] = size(data)

data = double(data);

% the size of the region should be at least 1 pixel and contain integer
% number of pixels
regionSize = max(round(regionSize),1)

% limit region size to the maximum number of pixels
regionSize = min(repmat(regionSize,1,3),[width, height, depth])

% initial centre on the distance from the edge
initStep = max(round(regionSize/2),1);

% determine centre coordinate for each dimension
% which should contain at least one point in each dimension
SX = initStep(1):regionSize(1):width;
SY = initStep(2):regionSize(2):height;
SZ = initStep(3):regionSize(3):depth;

% total initial number of clusters
K = length(SX) * length(SY) * length(SZ);

% array to store information about centres
% 1-color intensity
% 2-x coordinate
% 3-y coordinate
% 4-z coordinate
C = zeros(4,K);

CSize = size(C);

%% stop here - check

%initial coordinates of the centres on the equal-spaced grid
[XC,YC,ZC] = meshgrid(round(SX),round(SY),round(SZ));
C(2,:) = reshape(XC,1,K);
C(3,:) = reshape(YC,1,K);
C(4,:) = reshape(ZC,1,K);

%{
size(XC)
size(YC)
size(ZC)
size(C(2,:))
size(C(3,:))
size(C(4,:))
%}
%% stop here - check

% array to store original data
Xlab = zeros(width,height,depth,4);

% size(Xlab)
%% stop here - check

% assign coordinates of data
x = 1:width;
y = 1:height;
z = 1:depth;
[XL,YL,ZL] = meshgrid(y,x,z);

%size(XL)
%size(YL)
%size(ZL)

%% stop here2

Xlab(:,:,:,1) = data;
Xlab(:,:,:,2) = YL;
Xlab(:,:,:,3) = XL;
Xlab(:,:,:,4) = ZL;

% initial labels
l = ones(size(data))*-1;

% initial distance from the centres
d = inf(size(data));

% cluster centres moved in 3x3 neighborhood to the minimum intensity gradient
C = moveCentres3D(data, C, width, height, depth);

% threshold of distance between all centres
threshold = 0.5;

% initial residual error
E = inf;

% weighting factor for the spatial distance
w = 10;
sNorm = w/mean(regionSize);  

% values of C from previous iteration for comparison
C_old = C;

% number of pixels in each cluster (required for fast centres recomputation)
Ccount = zeros(1,K);
% mean value of colour intensity and coordinates of each cluster,
% which is required for fast centres recomputation.
Cmean = zeros(4,K);

% keep track of iterations
iteration = 0;

while E >= threshold
    for i = 1:K
        if (C(2,i)==0)
            continue;
        end
        % Calculate 2S subimage around each cluster
        minx = max(1,floor(C(2,i)-regionSize(1)));
        maxx = min(ceil(C(2,i)+regionSize(1)),width);
        miny = max(1,floor(C(3,i)-regionSize(2)));
        maxy = min(ceil(C(3,i)+regionSize(2)),height);
        minz = max(1,floor(C(4,i)-regionSize(3)));
        maxz = min(ceil(C(4,i)+regionSize(3)),depth);
        
        R2SX = minx:maxx;
        R2SY = miny:maxy;
        R2SZ = minz:maxz;
        dataLocal = Xlab(R2SX,R2SY,R2SZ,:);
        labelLocal = l(R2SX,R2SY,R2SZ);
        distLocal = d(R2SX,R2SY,R2SZ);
        dataLocal01 = dataLocal(:,:,:,1);
        dataLocal02 = dataLocal(:,:,:,2);
        dataLocal03 = dataLocal(:,:,:,3);
        dataLocal04 = dataLocal(:,:,:,4);
        dataLocal(:,:,:,1) = (dataLocal(:,:,:,1) - C(1,i));
        dataLocal(:,:,:,2) = (dataLocal(:,:,:,2) - C(2,i)).*sNorm;
        dataLocal(:,:,:,3) = (dataLocal(:,:,:,3) - C(3,i)).*sNorm;
        dataLocal(:,:,:,4) = (dataLocal(:,:,:,4) - C(4,i)).*sNorm;
        dataLocal = power(dataLocal,2);
        newDist = sum(dataLocal,4);
        selection = newDist < distLocal;
        updated = nnz(selection);
        if (updated)
            ls = labelLocal(selection);
            table = unique(ls);
            t1 = dataLocal01(selection);
            t2 = dataLocal02(selection);
            t3 = dataLocal03(selection);
            t4 = dataLocal04(selection);
            table(table==-1)=[];
            lt = length(table);
            sumQ = zeros(4,lt);
            sumQC = zeros(1,lt);
            for q = 1:lt
                qs = (ls == table(q));
                sumQC(q) = nnz(qs);
                sumQ(1,q) = sum(t1(qs));
                sumQ(2,q) = sum(t2(qs));
                sumQ(3,q) = sum(t3(qs));
                sumQ(4,q) = sum(t4(qs));
            end
            if (~isempty(table))
               Cmean(:,table) = ...
                   Cmean(:,table) - sumQ;
               Ccount(table) = Ccount(table) - sumQC;
            end
            Cmean(:,i) = Cmean(:,i) + [sum(t1); sum(t2); sum(t3); sum(t4)];
            Ccount(i) =  Ccount(i) + updated;
            
            distLocal(selection) = newDist(selection);
            labelLocal(selection) = i;
            d(R2SX,R2SY,R2SZ) = distLocal;
            l(R2SX,R2SY,R2SZ) = labelLocal;
        end
    end
    % update cluster centres via mean values
    
     selectionC = Ccount>0;
     C(:,~selectionC) = 0;
     C(:,selectionC) = (Cmean(:,selectionC)./repmat(Ccount(selectionC),4,1));
     
    % calculate residual error
    E = max(sqrt(sum((C_old(:,selectionC) - C(:,selectionC)).^2,1)));    
    C_old = C;
    iteration = iteration + 1;
    display(sprintf('Iteration:%d, residual error:%.2f, number of centres: %d',...
        iteration,E,length(selectionC)));
end

%%% Post-Processing for unlabeled pixels %%%

% remove labels that have no pixels
[newLabels,~,ic] = unique(l);
ind = find(newLabels == -1);
if (~isempty(ind))
    ic = ic - 1;
    ic(ic == ind) = -1;
    newLabels(ind) = [];
end
l(:) = ic;
K = length(newLabels);
C = C(:,newLabels);


% Create a mask for unlabeled pixels
I = (l == -1);

% Extract corresponding coordinates
XNL = XL(I);
YNL = YL(I);
ZNL = ZL(I);

% Loop through all selected points
for n = 1:length(XNL)
    xPoint = XNL(n);
    yPoint = YNL(n);
    zPoint = ZNL(n);
    point = data(yPoint,xPoint,zPoint); 
    
    % compute distance
    dc = point - C(1,:);
    dx = (xPoint - C(2,:)).*sNorm;
    dy = (yPoint - C(3,:)).*sNorm;
    dz = (zPoint - C(4,:)).*sNorm;
    newDist = sum(power([dc; dx; dy; dz],2),1);
    
    % select closest cluster as new label
    [~,ind] = min(newDist);
    l(xPoint,yPoint,zPoint) = ind; 
end
end
