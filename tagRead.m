function [header, tags] = tagRead(tagFile)
% This function reads in a (.tag) file.

% INPUT: name of .tag file

% OUTPUTS:
% header = header information of .tag file
% data = tagged pixel data

fileInfo = dir(tagFile);
% dir lists files and folders in the current folder (tagFile)
% e.g MyFolderInfo = dir('myfolder')
% MyFolderInfo = 5×1 struct array with fields:
%    name
%    folder
%    date
%    bytes
%    isdir
%    datenum

L = fileInfo.bytes;

% Read header information to string
% Open a file (tagFile)
fid = fopen(tagFile);

header = []; % empty array to store header info.
C = textscan(fid,'%s\n','delimiter',sprintf('\f'));
% textscan = read formatted data from text file or string
% %s\n = read series of characters, until find white space (%s), followed by new line (\n)
% sprintf = format data into string
% \f = form feed

for i = 1:length(C{1})-1 % we subtract 1 because there is an extra new line in the C array
    header = sprintf('%s', header,C{1}{i});
    % sprintf = format data into string
    % %s = character vector or string array.
    header = sprintf('%s\n', header);
    % %s = character vector or string array, followed by new line (\n).
end
header = sprintf('%s\f',header);
%s = character vector or string array,
% \f = form feed

% Go to last header position for start reading tags' pixels
fseek(fid, length(header), 'bof');
% fseek = Move to specified position in file
% fseek(fileID, offset, origin) sets the file position indicator offset bytes from origin in the specified file.
% 'bof' or -1 = beginning of file, in this case beginning of where header ends (last position)

% Read sizes of tags' data
headerV = strrep(header, ':', ' ');
% strrep = find and replace substring
% modifiedStr = strrep(origStr, oldSubstr, newSubstr)
% e.g. claim = 'This is a good example.';
% new_claim = strrep(claim, 'good', 'great')

C = textscan(headerV, '%s %s');
% textscan = read formatted data from text file or string
% %s %s = read series of characters, until find white space (two times in a row)

nameV = C{1,1};
valueV = C{1,2};

xx = str2double(cell2mat(valueV(strcmp(nameV, 'x'))));
yy = str2double(cell2mat(valueV(strcmp(nameV, 'y'))));
zz = str2double(cell2mat(valueV(strcmp(nameV, 'z'))));
% tf = strcmp(s1,s2) compares s1 and s2 and returns 1 (true) if 
% the two are identical and 0 (false) otherwise. Text is considered identical 
% if the size and content of each are the same. The return result tf is of data type logical.
% A = cell2mat(C) converts a cell array into an ordinary array. 
% The elements of the cell array must all contain the same data type, 
% and the resulting array is of that data type.

% Reading tags data
data = fread (fid, L - length(header), '*ubit8');
data = reshape(data,xx,yy,zz);
data = permute(data,[3 2 1]); % new rows = depth, new depth = rows
tags = permute(data,[2 3 1]); % new rows = columns, new columsn = depth, new depth = rows

% B = permute(A,order)
% Rearrange dimensions of N-D array
% B = permute(A,order) rearranges the dimensions of A so that 
% they are in the order specified by the vector order. 
% B has the same values of A but the order of the subscripts needed to
% access any particular element is rearranged as specified by order.
% All the elements of order must be unique, real, positive, integer values.
% For example:
% A = [1 2; 3 4]; permute(A,[2 1]) % switch rows with columns
% ans =
%     1     3
%     2     4

fclose(fid);

end