function [statisticsFeatures] = imageStatisticsFeatures(image)

% This function calculate the mean, median and standard deviation of a 2D image.
% imageStatisticsFeatures extracts the row vector containing mean, median and standard deviation from a 2D matrix.

% INPUT: 
% 2D image or part of an image in form of a 2D matrix.

% OUTPUT: 
% statisticsFeatures - row vector containing: 
% (a) mean, (b) median and (c) std (standard deviation) of 2D image.

      image = image(:); % transform to a single vector
      MEAN = mean(image); % calculate the mean of image
      MEDIAN = median(image); % calculate the median of image
      STD = std(image); % calculate the std. of image
      % return calculated characteristics
      statisticsFeatures = [MEAN MEDIAN STD];
end

