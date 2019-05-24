function data = readImageFile(filename)

% INPUT: file path to image data file (.img)
% OUTPUT: 3D array containing image data

    % read image metadata (data providing information about one or more aspects of the data)
    info = analyze75info(filename);
    % read data based on metadata
    data = analyze75read(info);
    % convert to same data type for all image files
    data = int16(data);
    % since Analyze 7.5 format uses radiological orientation (LAS),
    % flip the data for correct image display
    data = flip(data);
	
end