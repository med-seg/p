function [AllFeatures, AllCentres, AllTags] = ...
    patchFeatureExtractor(data, slicLabels, yPlus, tags)

% This function extracts features for fixed sized patches (25x25).
% This function extracts features and corresponding tags (annotated ground-truth) for ONE single .img (image data) file

% Inputs:
% data – 3D matrix loaded from .img (image data) file;
% slicLabels – corresponding SLIC labels;
% yPlus – estimation of the intensity probability density function estimated in pdfMap;
% tags – optional parameter: a 3D matrix of tags corresponding to the image data (if not provided .tag is assigned to -1 or unknown).

% Outputs: 
% AllFeatures - Nx46 matrix where each row consists of 46 features for a single (25x25) patch;
% AllCentres – Nx3 coordinate of the patch centre in the original image (x,y,z);
% AllTags – the label corresponding to the patch: 1 – if any of the pixels has tag > 0; -1 – if all tags are unknown; 0 – otherwise.


%{

%%%% begin: cropping of data and ground-truth (tags) to remove excess background %%%

tagCropped2 = [];

for ii = 1:50
 [BlobsBoundingBox,~] = cropping2(data,ii);
 tagSlice = tags(:,:,ii);
 tagSlice = imcrop(tagSlice, BlobsBoundingBox);
 tagSlice = imresize(tagSlice, [250 360]);
 tagCropped2(:,:,ii) = tagSlice;
end

tags = tagCropped2;


dataCropped = zeros(250,360,50);
for i = 1:50
    [thisBlobsBoundingBox,dataSlice] = cropping2(data,i);
    dataCropped(:,:,i) = dataSlice;
end
data = dataCropped;

%%%% end: cropping of data and ground-truth (tags) to remove excess background %%%

%}

if (nargin < 4)
    tags = -1*ones(size(data));
end

[width,height,slice_number] = size(data);

x = 3:3:width;
y = 3:3:height;
[XP,YP] = meshgrid(x,y);
XO = reshape(XP,numel(XP),1);
YO = reshape(YP,numel(YP),1);

AllFeatures = [];
AllCentres = [];
AllTags = [];

% loop through slices
for slc = 1:slice_number
    display(sprintf('patches features processing slice %d of %d: ',slc,slice_number));
    % read slice data
    slice = double(data(:,:,slc));
    tag_slice = tags(:,:,slc);
    
    [slice_modified, left, right, bottom, up ] = table_remove(slice);
	% slice_labels = slicLabels(left:right,bottom:up,slc);
    slice_labels = slicLabels(bottom:up,left:right,slc);
    tag_slice = tag_slice(bottom:up,left:right);
    
    % estimate probability distribution
    pdf_slice = map(slice_modified,yPlus);
    
    % calculate d_SIFT features for whole image slice
    [coords, d_patch] = vl_dsift(single(slice_modified),  'size', 6, 'geometry', [2,2,8]);
    xcoord  = coords(1,:);
    ycoord  = coords(2,:);
    
    %%% Remove centres outside of updated (cropped or cut) image slice %%%
    new_width = right - left + 1;
    new_height = up - bottom + 1;
    
    XR = XO - left + 1;
    YR = YO - bottom + 1;
    YR(XR<=12) = [];
    XR(XR<=12) = [];
    XR(YR<=12) = [];
    YR(YR<=12) = [];
    YR(XR>(right-left+1-12)) = [];
    XR(XR>(right-left+1-12)) = [];
    XR(YR>(up-bottom+1-12)) = [];
    YR(YR>(up-bottom+1-12)) = [];
    
    [nXR, nYR] = normalize_coordinate(XR, YR, [left, right], [bottom,up]);
    
    AllFeaturesSlice = zeros(length(XR),46);
    AllCentresSlice = zeros(length(XR),3);
    AllTagsSlice = zeros(length(XR),1);
    
    for i = 1:length(XR)
		% AllCentresSlice(i,:) = [XR(i), YR(i), slc];
          AllCentresSlice(i,:) = [YR(i) + bottom - 1, XR(i) + left - 1, slc];
        
        % initialise patch borders
        patchesX = max(1,XR(i)-12):min(XR(i)+12,new_width);
        patchesY = max(1,YR(i)-12):min(YR(i)+12,new_height);
        
        % read patch
        patch = slice_modified(patchesY,patchesX);
        pdf_patch = pdf_slice(patchesY,patchesX);
        slice_labels_patch = slice_labels(patchesY,patchesX);
        tag_patch = tag_slice(patchesY,patchesX);
        tag_patch = tag_patch(:);
        if (all(tag_patch<0))
            final_tag = -1;
        else
            if (any(tag_patch>0))
                final_tag = 1;
            else
                final_tag = 0;
            end
        end
        AllTagsSlice(i,:) = final_tag;
        
        % get region corresponding to superpixel of central pixel
        pixel_mask = slice_labels_patch==slice_labels_patch(13,13);
        pixel_region = patch(pixel_mask);
        pdf_region = pdf_patch(pixel_mask);
        
        
        idX = [1 find((xcoord <= (XR(i)+12))&(xcoord >= (XR(i)-12)) &...
            (ycoord <= (YR(i)+12))&(ycoord >= (YR(i)-12)))];
        c_current = coords(:,idX);
        [~,id] = min(sum(abs(c_current - repmat([XR(i);YR(i)],1,size(c_current,2)))));

        % descriptor closest to centre
        descriptor = d_patch(:,idX(id))';
        
        % patch statistics: mean, median, standard deviation
        patch_stat = imageStatisticsFeatures(patch);
        % probability statistics: mean, median, standard deviation
        pdf_stat = imageStatisticsFeatures(pdf_patch);
        
        
        region_stat = imageStatisticsFeatures(pixel_region);
        region_pdf_stat = imageStatisticsFeatures(pdf_region);
        
        AllFeaturesSlice(i,:) = [descriptor patch_stat pdf_stat region_stat region_pdf_stat nXR(i) nYR(i)];
    end
    
    AllFeatures = [AllFeatures; AllFeaturesSlice];
    AllCentres = [AllCentres; AllCentresSlice];
    AllTags = [AllTags; AllTagsSlice];
    
end

end

