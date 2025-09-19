function r_iau_enceladus = icrf2enceladus(r_icrf, r_enc, T)

    % Inputs
    % - r_icrf -> position of a body(ex. Sun) in ICRF
    % - T -> time in mjd2000
    %
    % Outputs
    % - r_iau_enceladus -> versor position of the same body(ex. Sun) in enceladus
    %                      body fixed reference frame(IAU Enceladus)

    % relative versor
    u_se = (r_icrf - r_enc) / norm(r_icrf - r_enc);

    % calculate enceladus orientation at specific time
    [RA, DEC, W] = enceladus_orientation(T);
    
    % From ICRF:
    % - z perpendicular to earth equator
    % - x direction earth vernal equinox
    % - y completes the orthornormal set) to Eceladus body frame 
    %
    % To Enceladus body fixed RF:
    % - z normal to Enceldus orbit around saturn (= enceladus equator)
    % - x points to the zero longitude direction
    % - y completes the orthornormal set)
    C = Rotz(W) * Rotx( (90 - DEC)) * Rotz((RA + 90)); 
    
    % Rotate from ICRF to Enceladus body fixed RF
    r_iau_enceladus = C * u_se';

end


%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Rotation Matrixes %%%
%%%%%%%%%%%%%%%%%%%%%%%%%
% These are the passive rotation matrix -> the rf frames are rotated
% This is why they have different sign from the classical active ones
% where the vectors are rotated and reference frame stays still

function R = Rotx(phi)
% INPUT -> phi [deg]

R = [1        0         0; ...
     0 cosd(phi)  sind(phi); ...
     0 -sind(phi)  cosd(phi)];
end

function R = Rotz(phi)
% INPUT -> phi [deg]

R = [cosd(phi)  sind(phi) 0; ...
     -sind(phi)  cosd(phi) 0; ...
       0           0      1];
end