function [r,v,kep] = EphSS_car_spice2(n, t, isSpice, spiceParam)
% Ephemerides of the solar system (SPICE‐enabled wrapper), without kernel
% loading

%%— NARGIN & DEFAULTS —————————————
if nargin < 3 || isempty(isSpice)
    isSpice = false;
end

if nargin < 4 || isempty(spiceParam)
    spiceParam.frame    = 'J2000';
    spiceParam.abcorr   = 'NONE';
    spiceParam.observer = '0';
end

%%— TWO‐BODY FALLBACK ——————————————
if ~isSpice
    mu  = 132724487690;  % Sun μ, km^3/s^2
    kep = uplanet(t,n);
    car = kep2car(kep,mu);
    r   = car(1:3);
    v   = car(4:6);
    return;
end 

% Compute NAIF ID string 
% spice_id = spice_ids(n);
idStr    = num2str(n);

% Compute state
et_time = UTC2ET(t); %conversion from utc to et(ephemeris time) (69s difference)
spice_time = et_time * 86400;
state      = cspice_spkezr(idStr, spice_time, spiceParam.frame,...
                           spiceParam.abcorr, spiceParam.observer);
state      = state';

r = state(1:3);
v = state(4:6);

% Compute Keplerians if requested
if nargout > 2
    if n < 10 || (mod(n,100)==0) || (mod(n,100)==99 && n<1000)
        mu = getAstroConstants('Sun','Mu');
    elseif n==301
        mu = getAstroConstants('Earth','Mu');
    elseif n>400 && n<499
        mu = getAstroConstants('Mars','Mu');
    elseif n>500 && n<599
        mu = getAstroConstants('Jupiter','Mu');
    elseif n>600 && n<699
        mu = getAstroConstants('Saturn','Mu');
    elseif n>700 && n<799
        mu = getAstroConstants('Uranus','Mu');
    elseif n>800 && n<899
        mu = getAstroConstants('Neptune','Mu');
    elseif n>900 && n<999
        mu = getAstroConstants('Pluto','Mu');
    else
        error('assign_mu:unknownID','No μ defined for ID %d.',n);
    end
    kep = car2kep(state,mu);
end
end
