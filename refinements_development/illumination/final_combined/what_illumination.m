function dataTable = what_illumination(dataTable, doEclipse, parallel)
% Calculates illumination conditions (terminator and eclipse) for given dataTable
% and adds eclipse and terminator columns to the existing dataTable
%
%   Inputs:
%       dataTable: Table with required columns:
%           - lat: latitude in radians
%           - lon: longitude in radians
%           - T: time in mjd2000
%       doEclipse (boolean): calculate or not eclipse conditions
%       parallel (logical): Flag for parallel processing
%
%   Output:
%       dataTable: Input table with added columns:
%           - terminator: 'night' or 'day'
%           - eclipse: 'Umbra', 'Penumbra', or 'no eclipse'

% Set defaults
if nargin < 3, parallel = false; end
if nargin < 2, doEclipse = false; end

% Validate input dataTable
if ~istable(dataTable)
    error('Input must be a table');
end
if ~all(ismember({'lat', 'lon', 'T'}, dataTable.Properties.VariableNames))
    error('dataTable must contain lat, lon, and T columns');
end

% Get the time value (assuming all T values are the same)
T = dataTable.T(1);

%% Initialize parameters for eclipse calculation
R_enc = 252.1;
spiceParam.frame = 'J2000';
spiceParam.abcorr = 'NONE';
spiceParam.observer = '602';

%% Calculate terminator (night/day)
dataTable.terminator = calculateTerminator(dataTable, T);

%% Calculate eclipse conditions 
if doEclipse
    % Only check eclipse for points on hemisphere facing Saturn (lon 0-90 and 270-360)
    eclipse_check_indices = find((dataTable.lon <= pi/2) | (dataTable.lon >= 3*pi/2));
    
    % Initialize all as 'no eclipse'
    dataTable.eclipse = repmat({'no eclipse'}, height(dataTable), 1);
    
    if ~isempty(eclipse_check_indices)
        % Create temporary table for eclipse calculation
        tempTable = dataTable(eclipse_check_indices, :);
        
        % Calculate eclipse conditions
        eclipse_conditions = calculate_lighting_conditions(tempTable, R_enc, spiceParam, parallel);
        
        % Map the string results back to our table
        for i = 1:length(eclipse_check_indices)
            idx = eclipse_check_indices(i);
            if strcmp(eclipse_conditions(i), 'Umbra')
                dataTable.eclipse{idx} = 'Umbra';
            elseif strcmp(eclipse_conditions(i), 'Penumbra')
                dataTable.eclipse{idx} = 'Penumbra';
            else
                dataTable.eclipse{idx} = 'no eclipse';
            end
        end
    end
else
    % No eclipse calculations 
    dataTable.eclipse = repmat({'no eclipse'}, height(dataTable), 1);
end

end

%% Helper function to calculate terminator (night/day)
function terminator = calculateTerminator(dataTable, T)
    % Get sub-solar point
    [lat_sun, lon_sun] = sunSubPointOnEnceladus(T);
    lon_sun = deg2rad(mod(lon_sun + 360, 360));
    
    % Convert to radians if needed
    lat_sun = deg2rad(lat_sun);
    
    % Calculate terminator for each point
    nPoints = height(dataTable);
    terminator = cell(nPoints, 1);
    
    for i = 1:nPoints
        lat_gt = dataTable.lat(i);
        lon_gt = dataTable.lon(i);
        
        % Calculate local hour angle
        lha = lon_gt - lon_sun;
        lha = mod(lha + pi, 2*pi) - pi; % Normalize to [-pi, pi]
        
        % Calculate terminator latitude at this longitude
        lat_term = atan(-cos(lha) / tan(lat_sun));
        
        % Determine if point is in night or day
        % Using the simple logic from original code
        if lat_sun >= 0
            % Sun to north: night is below terminator
            if lat_gt < lat_term
                terminator{i} = 'night';
            else
                terminator{i} = 'day';
            end
        else
            % Sun to south: night is above terminator
            if lat_gt > lat_term
                terminator{i} = 'night';
            else
                terminator{i} = 'day';
            end
        end
    end
end

