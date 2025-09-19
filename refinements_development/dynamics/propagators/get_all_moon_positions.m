function r_bodies_all_col = get_all_moon_positions(t_sec, t0_mjd, TB_IDs, spiceParam)
% Gets positions of specified bodies relative to a central body.
%
% INPUTS:
%   t_sec        (double) Time in seconds since the initial epoch.
%   t0_mjd       (double) The initial epoch as Modified Julian Date 2000.
%   TB_IDs       (vector) NAIF IDs of the target bodies.
%   spiceParam   (struct) SPICE parameters (.frame, .abcorr, .observer).

    % Convert time in seconds (from ode45) to an absolute MJD time
    current_mjd = t0_mjd + t_sec / 86400;
    
    num_TBs = length(TB_IDs);
    r_bodies_all_col = zeros(3, num_TBs); % Pre-allocate for column vectors

    % Define bodies that require special handling via the SSB
    special_handling_IDs = [5]; % 5 = Jupiter barycenter

    % Loop through each target body ID
    for i = 1:num_TBs
        target_id = TB_IDs(i);
        
        % This will hold the final position vector for the current body
        r_body_row = zeros(1, 3);
        
        % --- CONDITIONAL LOGIC ---
        if ismember(target_id, special_handling_IDs)
            % This body's position is not in the satellite kernel relative
            % to the central body. We must calculate it via the SSB.
            % R_target/central = R_target/SSB - R_central/SSB

            % Create a temporary SPICE param struct to override the observer to be the SSB
            spiceParamSSB = spiceParam;
            spiceParamSSB.observer = '0'; % Solar System Barycenter NAIF ID

            % 1. Get position of the target body relative to the SSB
            [r_target_wrt_ssb, ~] = EphSS_car_spice2(target_id, current_mjd, true, spiceParamSSB);
            
            % 2. Get position of the central body (e.g., Saturn) relative to the SSB
            idCentral = int32(str2double(spiceParam.observer));
            [r_central_wrt_ssb, ~] = EphSS_car_spice2(idCentral, current_mjd, true, spiceParamSSB);
            
            % 3. Perform the vector subtraction
            r_body_row = r_target_wrt_ssb - r_central_wrt_ssb;

        else
            % This is a local body (e.g., a Saturnian moon). We can query its
            % position directly relative to the central body (Saturn).
            % The original spiceParam.observer should be correct (e.g., '699').
            [r_body_row, ~] = EphSS_car_spice2(target_id, current_mjd, true, spiceParam);
        end
        
        % Transpose the 1x3 row vector to a 3x1 column vector and store it
        r_bodies_all_col(:, i) = r_body_row';
    end
end