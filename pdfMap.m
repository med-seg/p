function yPlus = pdfMap(N, Tagged, result_folder, isUpdate)

% This function computes the probability distributions of the pixels in image data.
% This estimation is performed via Gaussian kernel using a kernel density estimate for the positive and negative pixels.

% Inputs: 
% N – number of files to use as samples;
% Tagged – list of .img (image data) files to process
% result_folder – folder (path) where to save computation results for further computation 
% isUpdate - choose whether to recompute the results (1) load results from file (0)

% Output: yPlus


if (nargin < 4)
    isUpdate = true;
end

result_file = [result_folder '/PDF.mat'];
result_file = strrep(result_file,'//','/');

settings_file = [result_folder '/PDF_settings.mat'];
settings_file = strrep(settings_file,'//','/');

if (exist(result_file,'file') && isUpdate)
    delete(result_file);
else
    if(exist(result_file,'file'))
        load(result_file);
        return;
    end
end

if (exist(settings_file,'file'))
    load(settings_file);
else
    [min_pix, max_pix, ~, sigma_square] =...
        setupPDF(Tagged,settings_file,isUpdate);
end

samples = randperm(length(Tagged),min(N,length(Tagged)));

sigma_sqd = -1/2*sigma_square;

H = double(min(0,min_pix):max(max_pix,512))';
range_length = length(H);

% estimate positive pixels
FP = zeros(range_length,1);
% estimate negative pixels
FM = zeros(range_length,1);

% count of positive class pixels
n = 0;
% count of negative class pixels
m = 0;

% number of pixels to proceess simultaneously to avoid memory overflow
part = 100000;

for i = 1:length(samples)
    % path to the original file with data
    datafile = [Tagged(samples(i)).path Tagged(samples(i)).name '.img'];    
    tagfile = [Tagged(samples(i)).path Tagged(samples(i)).name '.tag'];
    
    display(sprintf('PDF processing file %d of %d: ',i,length(samples)));
    
    data = readImgFile(datafile); % read .img (image datat) file into a 3D array
    [~, tags] = tagRead(tagfile); % read .tag (ground-truth of manual annotation) file into a 3D array

    
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


    x_plus = tags > 0;
    x_minus = ~x_plus;
    
    XP = cast(data(x_plus),'double');    
    XM = cast(data(x_minus),'double');
    psize = length(XP);
    msize = length(XM);
    n = n + psize;
    m = m + msize;
    part_numberP = ceil(psize/part);
    part_numberM = ceil(msize/part);

    for j = 1:part_numberP
        pos = 1+(j-1)*part:min(psize,j*part);
        FP = FP + sum(exp((pdist2(H, XP(pos)).^2)*sigma_sqd),2);
    end
    
    for j = 1:part_numberM      
        display(sprintf('Part %d of %d: ',j,part_numberM));
        pos = 1+(j-1)*part:min(msize,j*part);       
        FM = FM + sum(exp((pdist2(H, XM(pos)).^2)*sigma_sqd),2);
    end
end

FP = FP/n;
FM = FM/m;

yPlus = FP./(FP+FM);
yPlus(isnan(yPlus)) = 0;
save(result_file,'yPlus','-mat');
end

