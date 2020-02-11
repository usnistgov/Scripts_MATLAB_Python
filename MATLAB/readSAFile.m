function [data5,dataType] = readSAFile(filename,delimiter,TIME_COL,TEMP_COL,RMS_COL)
% [data5,dataType] = readFaroDataFormat4(filename,delimiter,TIME_COL,TEMP_COL)
% This function requires atleast the 'filename' and 'delimiter' and it will
% return the coordinates in the file. TIME and TEMPERATURE data are
% optional
% It returns the data, the time, temperature and the dataType, which could
% be spherical, cartesian, or cylindrical
%
% Return data will have x,y,z (or r,th,phi) and the time,temperature and
% rms values

% This is standard start row, column and the # of columns.
start_row = 5;
start_col = 2;
num_of_columns = 3;

if (nargin<3)
    TIME_COL = 'NONE';
    TEMP_COL = 'NONE';
    RMS_COL = 'NONE';
end
if (nargin<4)
    TEMP_COL = 'NONE';
    RMS_COL = 'NONE';
end
if (nargin<5)
    RMS_COL = 'NONE';
end


data = importdata(filename,delimiter);
tline1 = char(data{end,:});
dateCol = findDateCol(tline1,delimiter);
tempCol = findTempCol(tline1,delimiter,'Weather:T=');
rmsCol = findRMSCol(tline1,delimiter,'RMS Error');
%tline2 = char(data{end,:});
%fprintf('Date column = %d\n',dateCol);
FLAG = 1;

%Check the first 20 lines if there's any information about the data type.
for jj = 1:min(20,length(data))
    tline = char(data{jj,:});
    
    if  strcmp(tline(1:19), '// AXES = SPHERICAL' )
        dataType = 'Spherical';
        break;
    elseif strcmp(tline(1:19), '// AXES = CARTESIAN')
        dataType = 'Cartesian';
        break;
    elseif strcmp(tline(1:19), '// AXES = CYLINDRIC')
        dataType = 'Cylindrical';
        break;
    else
        dataType = 'Unknown';
    end
end

for jj=start_row+1:length(data)
    tline = char(data{jj,:});
    
    
    
    idx =jj-start_row;
    try
    data3(idx,:) = regexp(tline, delimiter, 'split');
    catch
        error('It is likely that the exported file does not have the same number of columns');
    end
    data4 = str2num(char(data3{idx,start_col:start_col+num_of_columns-1}));
    if( strcmp(TIME_COL,'NONE')~=1 && dateCol >0 )
        dateNumber = datenum(datestr(data3{idx,dateCol}));
    else
        dateNumber = idx;
    end
    
    if( strcmp(TEMP_COL,'NONE')~=1 && tempCol >0 )
        weatherData1 = char(data3{idx,tempCol});
        aa = strsplit(weatherData1,'='); bb = aa{2};
        tempData = str2double(bb(1:end-1));
        tempData = (tempData - 32)/1.8;
    else
        %disp('No temp col');
        tempData = 0;
    end
    
    if( strcmp(RMS_COL,'NONE')~=1 && rmsCol >0 )
        rmsData1 = char(data3{idx,rmsCol});
        aa2 = strsplit(rmsData1,' '); bb2 = aa2{4};
        rmsData = str2double(bb2);
    else
        %disp('No temp col');
        rmsData = 0;
    end
    
    
    data5(idx,:) = [data4' dateNumber tempData rmsData];
end


function dateCol = findDateCol(tline,delimiter)
%tline = char(data{end,:});
allCols = regexp(tline, delimiter, 'split');
for i = 1:length(allCols)
    slashes = length(strfind( (allCols{i}),'/'));
    colons = length(strfind( (allCols{i}),':'));
    if slashes == 2 && colons == 2
        dateCol = i; %disp('Time found')
        break;
    else
        dateCol = -1;
    end
end

function tempCol = findTempCol(tline,delimiter,findString)
%tline = char(data{end,:});
allCols = regexp(tline, delimiter, 'split');
for i = 1:length(allCols)
    tempData = length(strfind( (strtrim(char(allCols{i}))),'Weather:T='));
    %char(allCols{i})
    if tempData == 1
        tempCol = i; %disp('Temp found')
        break;
    else
        tempCol = -1;
    end
end

function tempCol = findRMSCol(tline,delimiter,findString)
%tline = char(data{end,:});
allCols = regexp(tline, delimiter, 'split');
for i = 1:length(allCols)
    tempData = length(strfind( (strtrim(char(allCols{i}))),findString));
    %char(allCols{i})
    if tempData == 1
        tempCol = i; %disp('Temp found')
        break;
    else
        tempCol = -1;
    end
end