function jd2000_tdb = UTC2ET(jd2000_utc)
%  Manually converts J2000 UTC, the time used in GIGANTES to J2000 TDB
%  (ET), the time used by SPICE. 
% 
% On Enceladus they differ averagely of 70s.
%
%  Based on the SPICE DELTET routine and naif00012.tls leapseconds kernel
%  values.
%   
% INPUT:
%   jd2000_utc[1]  in days
%
% OUTPUT:
%   jd2000_tdb[1]  in days


    %% Part 1: Hard-coded constants from the LSK file
    DELTA_T_A = 32.184;
    K         = 1.657e-3;
    EB        = 1.671e-2;
    M0        = 6.239996;
    M1        = 1.99096871e-7;

    % DELTA_AT table: [leap_seconds, MJD2000_UTC_of_change]
    % The dates from the LSK have been pre-converted to MJD2000 UTC
    % for direct comparison.
    leap_second_table_mjd = [
        10, -10224.0; % 1972-JAN-01 00:00:00 UTC -> -10223.5 MJD -> -10224.0 MJD2000
        11, -10042.0; % 1972-JUL-01
        12,  -9858.0; % 1973-JAN-01
        13,  -9493.0; % 1974-JAN-01
        14,  -9128.0; % 1975-JAN-01
        15,  -8763.0; % 1976-JAN-01
        16,  -8397.0; % 1977-JAN-01
        17,  -8032.0; % 1978-JAN-01
        18,  -7667.0; % 1979-JAN-01
        19,  -7302.0; % 1980-JAN-01
        20,  -6755.0; % 1981-JUL-01
        21,  -6390.0; % 1982-JUL-01
        22,  -6025.0; % 1983-JUL-01
        23,  -5294.0; % 1985-JUL-01
        24,  -4380.0; % 1988-JAN-01
        25,  -3649.0; % 1990-JAN-01
        26,  -3284.0; % 1991-JAN-01
        27,  -2737.0; % 1992-JUL-01
        28,  -2372.0; % 1993-JUL-01
        29,  -2007.0; % 1994-JUL-01
        30,  -1458.0; % 1996-JAN-01
        31,  -1093.0; % 1997-JUL-01
        32,   -362.0; % 1999-JAN-01
        33,   2193.0; % 2006-JAN-01
        34,   3289.0; % 2009-JAN-01
        35,   4566.0; % 2012-JUL-01
        36,   5661.0; % 2015-JUL-01
        37,   6211.0; % 2017-JAN-01
    ];
    
    %% Part 2: Calculate DELTA_ET using an iterative approach
    
    % Convert input MJD UTC to seconds past J2000 UTC
    utc_seconds_past_j2000 = jd2000_utc * 86400.0;
    
    % Initial guess for ET seconds past J2000
    et_seconds_approx = utc_seconds_past_j2000;
    
    % Iterate a few times to converge (2 is usually enough)
    for i = 1:3
        % Find the correct number of leap seconds (DELTA_AT) for our date
        idx = find(leap_second_table_mjd(:, 2) <= jd2000_utc, 1, 'last');
        if isempty(idx)
            % Should not happen for modern dates
            DELTA_AT = 10;
        else
            DELTA_AT = leap_second_table_mjd(idx, 1);
        end
        
        % Calculate Mean Anomaly (M)
        t = et_seconds_approx;
        M = M0 + M1 * t;
        
        % Calculate Eccentric Anomaly (E) via iteration
        E = M;
        for j = 1:5
            E = M + EB * sin(E);
        end
        
        % Calculate ET - TAI
        ET_minus_TAI = DELTA_T_A + K * sin(E);
        
        % Calculate DELTA_ET
        DELTA_ET = ET_minus_TAI + DELTA_AT;
        
        % Refine our approximation of ET
        et_seconds_approx = utc_seconds_past_j2000 + DELTA_ET;
    end
    
    % Convert the final ET in seconds back to JD2000 TDB
    jd2000_tdb = et_seconds_approx / 86400.0;
end