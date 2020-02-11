function [fileNames] = getSortedFiles(folder1, fileType, sortType)
% Get desired files by sortType (datetime, size or name)
files1 = dir(fullfile(folder1,fileType));

if strcmp(sortType,'datetime')==1
    parameters1 = [files1(:).datenum]';
    [~,idx1] = sort(parameters1);
    fileNames = {files1(idx1).name};
    
elseif strcmp(sortType,'size')==1
    parameters1 = [files1(:).bytes]';
    [~,idx1] = sort(parameters1);
    fileNames = {files1(idx1).name};
    
elseif strcmp(sortType,'name')==1
    parameters1 = {files1.name};
    [~,idx1] = sort(parameters1);
    fileNames = {files1(idx1).name};
end