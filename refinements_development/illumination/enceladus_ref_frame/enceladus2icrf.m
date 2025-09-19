function r_icrf = enceladus2icrf(u_ebf, r_enc, R_enc, T)

    % Inputs
    % - u_ebf -> row versor in EBF - Enceladus Body Fixed reference frame
    % - r_enc -> enceladus position in ICRF
    % - R_enc -> spherical radius Enceladus
    % - T -> time in mjd2000
    %
    % Outputs
    % - r_icr -> position in ICRF
    
    % check if the u_ebf is a versor
    tol = 1e-6;
    if abs(norm(u_ebf) - 1) < tol
        u_ebf = (u_ebf'/norm(u_ebf))*R_enc;
    else % project the versor on the Enceladus surface(assumed spherical)
        u_ebf = u_ebf'*R_enc;
    end

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
    
    % position of the point from the center of enceladus in ICRF
    r_fromthecenter_icrf = C' * u_ebf;

    % total position in ICRF 
    r_icrf = r_enc +  r_fromthecenter_icrf';
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