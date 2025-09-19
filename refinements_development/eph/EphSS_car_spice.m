function [r,v,kep] = EphSS_car_spice(n, t, isSpice, kernels, spiceParam)
% EphSS_car_spice.m – Ephemerides of the solar system (SPICE‐enabled wrapper)
% Two‐phase design with advanced per-session kernel management.

%%— NARGIN & DEFAULTS —————————————
if nargin < 3 || isempty(isSpice)
    isSpice = false;
end
if nargin < 4 || isempty(kernels)
    kernels = {'sat441.bsp','naif0012.tls'};
    if ismember(n, [1,2,4,5,7,8,9])
        isSpice = false;
    end
end
if nargin < 5 || isempty(spiceParam)
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

%%— SPICE BRANCH WITH PERSISTENT INITIALIZATION ———————————
% Store the state of the last successful initialization
persistent lastKernels lastSpiceParam initialized;

% --- MODIFIED LOGIC ---
% Determine if re-initialization is needed.
% It's needed if:
% 1. It's the first run (initialized is empty).
% 2. The kernel list has changed.
% 3. The SPICE parameters have changed.
needs_init = isempty(initialized) || ...
             ~isequal(kernels, lastKernels) || ...
             ~isequal(spiceParam, lastSpiceParam);

if needs_init
    fprintf('SPICE RE-INITIALIZING: Kernels or parameters changed.\n');

    % Clear all previously loaded kernels to ensure a clean slate
    cspice_kclear; 
    
    % Resolve paths and load all kernels
    kernelsToLoad = {};
    for i = 1:numel(kernels)
        kernelFile = string(kernels{i});
        kernelPath = which(kernelFile);
        assert(~isempty(kernelPath), 'Cannot find kernel file "%s".', kernelFile);
        kernelsToLoad{end+1} = kernelPath;
    end
    
    % Load all resolved kernels
    cspice_furnsh(kernelsToLoad);

    % Store the current kernels and parameters for the next call
    lastKernels    = kernels;
    lastSpiceParam = spiceParam;
    initialized    = true;
end


% Compute NAIF ID string each call
spice_id = spice_ids(n);
idStr    = num2str(spice_id);

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
