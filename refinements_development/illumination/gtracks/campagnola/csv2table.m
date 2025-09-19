clear
clc

csvPath = which('mergedTables.csv');
if isempty(csvPath)
    error('Cannot locate mergedTables.csv – please add its folder to the MATLAB path.');
end
filename = csvPath;

% 1) Read all lines
fid = fopen(filename,'r');
if fid<0
    error('Cannot open file: %s', filename);
end
lines = textscan(fid,'%s','Delimiter','\n','Whitespace','');
fclose(fid);
lines = lines{1};
if numel(lines)<2
    error('File has no data rows');
end

% 2) Split header, count columns
headerFields = strsplit(lines{1},',');
nCols = numel(headerFields);

% 3) Parse each data line, check field count
nData = numel(lines)-1;
dataCells = cell(nData, nCols);
for i = 2:numel(lines)
    rowIdx = i-1;
    fld = strsplit(lines{i},',');
    if numel(fld)~=nCols
        error('Line %d: found %d fields, expected %d', i, numel(fld), nCols);
    end
    dataCells(rowIdx,:) = fld;
end
fprintf('All %d data rows have %d fields → formatting OK\n', nData, nCols);

% 4) Make valid MATLAB names for the variables
varNames = matlab.lang.makeValidName(headerFields);

% 5) Build table (everything still as strings)
T = cell2table(dataCells,'VariableNames',varNames);

% 6) Convert these columns to numeric
numericCols = { ...
  'ALT','v_inf','v', ...
  'LAT','LON','SEP','PHS','TRA','INC', ...
  'PER','n','m','dTID','TID', ...
  'dDeltaV','DeltaV','dTOF','TOF','SJB', ...
  'kappa_minus','kappa_plus','alpha_minus','alpha_plus', ...
  'RP','RA','ECLIP','ETRA','SIIILON' ...
};
for ii = 1:numel(numericCols)
    col = numericCols{ii};
    T.(col) = str2double(T.(col));
end
