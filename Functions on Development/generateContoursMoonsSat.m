function [rascCONT, rpscCONT, alpha] = generateContoursMoonsSat(pl, vInf, idpl)

% only works for elliptical orbits
% Saturn is the main body 

if idpl == 1
    mu = 132724487690;
    [~, ~, rPL] = planetConstants(pl);
elseif idpl == 5
    mu  = planetConstants(idpl);
    rPL = jupMoonsConstants(pl);
elseif idpl == 6
    mu  = planetConstants(idpl);
    rPL = satMoonsConstants(pl);
elseif idpl == 7
    mu  = planetConstants(idpl);
    rPL = uranusMoonsConstants(pl);
end
alpha = deg2rad(linspace(0,180));

vPL = sqrt(mu/rPL);

rpscCONT = zeros(length(alpha),1);
rascCONT = zeros(length(alpha),1);

for indi = 1:length(alpha)
    [rascCONT(indi,:), rpscCONT(indi,:)] = SCorbitMoonsSat(alpha(indi), vInf, vPL, rPL);
end

% --> eliminate hyperbolic orbits
idxs           = find(rascCONT < 0);
rascCONT(idxs) = [];
rpscCONT(idxs) = [];
alpha(idxs)    = [];

% --> eliminate retrograde orbits
[~, row] = min(rpscCONT);
rascCONT(row+1:end) = [];
rpscCONT(row+1:end) = [];
alpha(row+1:end)    = [];

end
