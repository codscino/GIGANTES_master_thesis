function lighting_conditions = calculate_lighting_conditions(dataTable, R_enc, spiceParam, parallel, kernels)
%
%  This function determines if points of a groundtrack on Enceladus are
%  in sunlight, umbra, or penumbra cast by Saturn. 
% 
%  It calculates only once at the flyby pericentre time(half time):
%   - apparent size of the Sun
%   - beta -> sun elevation angle from enceladus 
%   - apparent angular size of spheroid Saturn(in function of beta)
%   - sun ICRF position
%  
%  Then it iterates through each data point to determine the lighting
%
%  Author : Claudio Ferrara
%
%  Inputs:
%    dataTable : A MATLAB table with required fields: 'lat', 'lon', 'T'
%                       'lat' and 'lon' are in [rad], 'T' is J2000 ET time.
%
%  Outputs:
%    lighting_conditions : A cell array of strings ('Lit', 'Umbra', 'Penumbra')
%                          corresponding to each row of the input table.

    if nargin < 3
        spiceParam.frame    = 'J2000';
        spiceParam.abcorr   = 'NONE';
        spiceParam.observer = '0';
    end
    if nargin < 4
        parallel = false;
    end

    if nargin < 5
        kernels = {'sat441.bsp', 'naif0012.tls'};
    end

    if parallel
        loadSpiceKernels(kernels, parallel) % load kernels for every parallel worker
    end


    %% Step 1: Calculate "Constant" Values Once
    % Use the middle row of the table(pericentre) as the representative point in time
    mid_idx = round(height(dataTable) / 2);
    T_mid = dataTable.T(mid_idx);
    
    % Get body positions at this representative time
    r_enc_mid = EphSS_car_spice2(602, T_mid, true, spiceParam);
    r_sat_mid = EphSS_car_spice2(699, T_mid, true, spiceParam);
    r_sun_mid = EphSS_car_spice2(10,  T_mid, true, spiceParam);

    % Calculate the apparent angular DIAMETER of spherical Sun from
    % Enceladus center
    alpha_sun = getAngSize(r_enc_mid, r_sun_mid, 10);

    %%% Calculate the apparent angular DIAMETER of  spheroid Saturn from Enceladus surface %%%
    % Distance between Enceladus and Saturn centers
    d = norm(r_enc_mid - r_sat_mid);
    % Averaged distance from Enceladus surface(facing Saturn) to Saturn center
    d_mean = emisphere_distances(d, R_enc); 
    % Enceladus declination angle between Sun and Saturn
    beta = saturn_sun_beta_angle(r_sun_mid, r_sat_mid, r_enc_mid, T_mid);
    % true angular size of spheroid Saturn 
    alpha_sat = saturnEclipse(beta, d_mean);
   
   
    %% Step 2: Iterate Through Each Point and Determine Lighting
    num_points = height(dataTable);
    lighting_conditions = cell(num_points, 1);

    if parallel
        parfor i = 1:num_points
            T_i   = dataTable.T(i);
            lat_i = dataTable.lat(i);
            lon_i = dataTable.lon(i);

            r_enc_i = EphSS_car_spice2(602, T_i, true, spiceParam);
            r_sat_i = EphSS_car_spice2(699, T_i, true, spiceParam);

            sep_rad = angle_sep(lat_i, lon_i, T_i, r_enc_i, r_sat_i, r_sun_mid, R_enc);

            sum_of_radii = alpha_sat + alpha_sun;
            diff_of_radii = abs(alpha_sat - alpha_sun);

            if sep_rad <= diff_of_radii
                lighting_conditions{i} = 'Umbra';
            elseif sep_rad > diff_of_radii && sep_rad < sum_of_radii
                lighting_conditions{i} = 'Penumbra';
            else
                lighting_conditions{i} = 'Lit';
            end
        end
    else
        for i = 1:num_points
            T_i   = dataTable.T(i);
            lat_i = dataTable.lat(i);
            lon_i = dataTable.lon(i);

            r_enc_i = EphSS_car_spice2(602, T_i, true, spiceParam);
            r_sat_i = EphSS_car_spice2(699, T_i, true, spiceParam);

            sep_rad = angle_sep(lat_i, lon_i, T_i, r_enc_i, r_sat_i, r_sun_mid, R_enc);

            sum_of_radii = alpha_sat + alpha_sun;
            diff_of_radii = abs(alpha_sat - alpha_sun);

            if sep_rad <= diff_of_radii
                lighting_conditions{i} = 'Umbra';
            elseif sep_rad > diff_of_radii && sep_rad < sum_of_radii
                lighting_conditions{i} = 'Penumbra';
            else
                lighting_conditions{i} = 'Lit';
            end
        end
    end
end