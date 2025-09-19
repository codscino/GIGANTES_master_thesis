function spice_id = spice_ids(uplanet_id)
% Function that converts from uplanet numbering to Spice NAIF numbering.
%
% uplanet numebring:
% 1:   Mercury
% 2:   Venus
% 3:   Earth
% 4:   Mars
% 5:   Jupiter
% 6:   Saturn
% 7:   Uranus
% 8:   Neptune
% 9:   Pluto
% 10:  Sun
%
% Spice NAIF numbering:
% from 1 to 9 spice instead of giving the planet position gives the planet
% baricentre position. To get planet position add a 99:
% 199:   Mercury
% 299:   Venus
% 399:   Earth
% 499:   Mars
% 599:   Jupiter
% 699:   Saturn
% 799:   Uranus
% 899:   Neptune
% 999:   Pluto

if uplanet_id >= 1 && uplanet_id <= 9
    spice_id = 100*uplanet_id + 99;
else 
    spice_id = uplanet_id;
end

% for Sun 10 works is the same between uplnaet and spice,
% for the solar system baricentre SSB, spice uses 0

% the moon of Saturn follow the Spice NAIF numeration:
% 
% MIMAS (601) w.r.t. SATURN BARYCENTER (6)
% ENCELADUS (602) w.r.t. SATURN BARYCENTER (6)
% TETHYS (603) w.r.t. SATURN BARYCENTER (6)
% DIONE (604) w.r.t. SATURN BARYCENTER (6)
% RHEA (605) w.r.t. SATURN BARYCENTER (6)
% TITAN (606) w.r.t. SATURN BARYCENTER (6)
% HYPERION (607) w.r.t. SATURN BARYCENTER (6)
% IAPETUS (608) w.r.t. SATURN BARYCENTER (6)
% PHOEBE (609) w.r.t. SATURN BARYCENTER (6)
% HELENE (612) w.r.t. SATURN BARYCENTER (6)
% TELESTO (613) w.r.t. SATURN BARYCENTER (6)
% CALYPSO (614) w.r.t. SATURN BARYCENTER (6)
% METHONE (632) w.r.t. SATURN BARYCENTER (6)
% POLYDEUCES (634) w.r.t. SATURN BARYCENTER (6)


end