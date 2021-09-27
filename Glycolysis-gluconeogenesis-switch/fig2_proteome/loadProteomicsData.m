%% Import data from csv file.

%% Initialize variables.
filename = '2015_10_29_lag_phase_ratios.csv';
delimiter = ',';
startRow = 2;

%% Format string for each line of text:
%   column1: double (%f)
%	column2: text (%s)
%   column3: double (%f)
%	column4: text (%s)
%   column5: double (%f)
%	column6: double (%f)
%   column7: text (%s)
%	column8: text (%s)
formatSpec = '%f%s%f%s%f%f%s%s%[^\n\r]';

%% Open the text file.
fileID = fopen(filename,'r');

%% Read columns of data according to format string.
dataArray = textscan(fileID, formatSpec, 'Delimiter', delimiter, 'EmptyValue' ,NaN,'HeaderLines' ,startRow-1, 'ReturnOnError', false);

%% Close the text file.
fclose(fileID);

%% Allocate imported array to column variable names
VarName1 = dataArray{:, 1};
protein = dataArray{:, 2};
charge = dataArray{:, 3};
seq = dataArray{:, 4};
svmPred = dataArray{:, 5};
ratio = dataArray{:, 6};
condition = dataArray{:, 7};
gene = dataArray{:, 8};

%% Clear temporary variables
clearvars filename delimiter startRow formatSpec fileID dataArray ans;