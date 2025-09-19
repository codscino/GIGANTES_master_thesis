function r0_sc_offset = offsetR0_sc(r0_sc, pars)

%   Inputs:
%     r0_sc  - 3×1, initial saturn centric spacecraft position in km
%     pars   - struct with fields:
%                pars.INPUTS.Flyby.min_h  (double) minimum flyby altitude in km
%
%   Outputs:
%     r0_sc_offset  - 3×1, offset spacecraft position in km, saturn
%     centric
 
    % Apply flyby‐minimum‐altitude offset to a spacecraft position.

    % unit vector from the central body toward the spacecraft
    r_enc_unit = r0_sc ./ norm(r0_sc);

    % build an in‐plane, perpendicular direction (rotated by +90° in the XY‐plane)
    r_offset_direction = [-r_enc_unit(2); r_enc_unit(1); 0];

    % flyby pericentre altitude, 252km is enceladus radius
    rp_flyby = pars.INPUTS.Flyby.min_h + 252;

    % offset vector and apply
    r_offset = r_offset_direction * rp_flyby;
    r0_sc_offset = r0_sc + r_offset';
end