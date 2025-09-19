function spiceCheck(time, id, kernels)

if nargin < 3
    % use the saturn kernel
    SPKfile = 'sat441.bsp';
    leapsecondsKernel = 'naif0012.tls';
    kernels = {SPKfile, leapsecondsKernel};
end

spice_id = spice_ids(id);

SPKfile = string(kernels(1));
leapsecondsKernel = string(kernels(2));
spkname =  which(SPKfile);
LSK   = which(leapsecondsKernel);
cspice_furnsh({spkname,LSK});

kernelRange  = cspice_spkcov(spkname, spice_id, 1);

if isempty(kernelRange)
    error(['The requested body (ID %d) is not available in the kernel "%s". ' ...
           'Try downloading and selecting another SPK from:\n' ...
           'https://naif.jpl.nasa.gov/pub/naif/generic_kernels/spk/'], ...
           id, SPKfile);
end

spice_startDate = kernelRange(1)/86400; % seconds -> days
spice_endDate = kernelRange(2)/86400; % seconds -> days

% extract dates from the time vector
if length(time) > 1
    startDate = time(1);
    endDate = time(2);
else
    startDate = time(1);
end
    

if length(time) < 2
     % ---------------------------------------------------
     % single date requested: t = startDate
     % ---------------------------------------------------
    if startDate < spice_startDate || startDate > spice_endDate
        % out of bounds
        warning('Date out of bounds for available SPK file.')
        fprintf('Coverage for object %d\n', id )
        timstr = cspice_timout(kernelRange','YYYY MON DD HR:MN:SC (TDB) ::TDB');
        fprintf('Start: %s\n', timstr(1,:));
        fprintf('Stop: %s\n', timstr(2,:));
        error('Please insert a date between the times above or use another SPK file.');
    end

else
    % ---------------------------------------------------
    % date span requested: t = startDate:endDate
    % ---------------------------------------------------
    if startDate < spice_startDate
        timstr = cspice_timout(kernelRange', 'YYYY MON DD HR:MN:SC (TDB) ::TDB');
        fprintf('Earliest available for object %d: %s\n', id, timstr(1,:));
        error('Please choose a startDate ≥ the earliest coverage time.');
    elseif endDate > spice_endDate
        timstr = cspice_timout(kernelRange', 'YYYY MON DD HR:MN:SC (TDB) ::TDB');
        fprintf('Latest available for object %d:   %s\n', id, timstr(2,:));
        error('Please choose an endDate ≤ the latest coverage time.');
    elseif startDate > endDate
        error('startDate must be ≤ endDate.');
    end

end


%cspice_kclear %clear kernel
% cspice_unload({ spkname, leapsecondsKernel }); %unloads kernels

end