**IMPORTANT**: You are encouraged to read this entire file and make necessary modifications to your code before program execution.

MATLAB code for [description of paper]

by Hykoush Asaturyan, Barbara Villarini, E. Louise Thomas, Jimmy. D. Bell and Julie Fitzpatrick.

Please see [link to pdf file] for online access to PDF file.

This program was built using MATLAB R2017b and MATLAB R2018a, and employs the Deep Learning Toolbox, Statistics and Machine Learning Toolbox, Image Processing Toolbox and Computer Vision Toolbox.

**1) TRAINING STAGE:**

To initiate the training stage, modify and run "runFunction.m".

**IMPORTANT TO NOTE:**  Please modify the following variables below to suit image data file type, ground-truth (annotation) file type and options for training the deep learning model. Replace 'define' with appropriate value.

imWidth = 'define';

imHeight = 'define';

maxEpoch = 'define';  %300, 320, 420

batchSize = 'define'; % 1, 10

**IMPORTANT TO NOTE:** "cropping2.m" is a function that aims to remove background (black) and crop to contain abdomen scan. If necessary, please modify or use an appropriate cropping function for your image data and apply the same cropping to the ground-truth labels.

**2) TESTING STAGE:**

To initiate the training stage, modify and run one of the following files:

**generate_predictions.m** -  generate binary segmentation mapping using trained deep learning model. You will need to modify this file accordingly.

**generate_limit_range_predictions.m** -  generate binary segmentation mapping using trained deep learning model across probability thresholds [0.05,0.95]. You will need to modify this file accordingly.

**IMPORTANT TO NOTE:** Please modify the following variables below to suit image data file type and ground-truth (annotation) file type. Replace 'define' with appropriate value.

imWidth = 'define';

imHeight = 'define';

**IMPORTANT TO NOTE:** For digital contrast enhancement (DCE), it is important to train an artifical neural network to predict the DCE function parameters of gain and cut-off value. The purpose of the DCE is to increase intensity variation between pancreas and non-pancreas pixels in the test image (as a volumetric 3D array) on a slice-by-slice basis.

There are two possible options to use this stage.

**(a) artifical neural network (6 inputs; 2 hidden layers, each layer has 20 nodes; 2 outputs).**

**6 inputs include the mean intensity of:**
1) 3D array containing entire image data;
2) middle slice (2D image in 3D array); 
3) slice at the first quartile in 3D array;
4) slice at the third-quartile in 3D array;
5) 3D array containing first quartile of image data (e.g. for a 3D array containing 80 slices, the first quartile would contain the first 20 slices).
6) 3D array containing last quartile of image data (e.g. for a 3D array containing 80 slices, the third quartile would contain the last 20 slices).

**2 outputs:**
1) desired gain
2) desired cut-off value

**(b) artifical neural network for predicitng gain AND  linear regression function for predicting cut-off value based on given gain.**

1) artifical neural network (6 inputs from (a); 2 hidden layers, each layer has 20 nodes; 1 output for gain).
2) linear regression model (input: desired gain; output: desired cut-off value).




