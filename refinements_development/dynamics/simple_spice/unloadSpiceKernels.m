function unloadSpiceKernels(kernels)
%   Unload SPICE kernels using WHICH + cspice_unload.
%
%   Example:
%       kernels = {'sat441.bsp', 'naif0012.tls', 'pck00010.tpc'};
%       unloadSpiceKernels(kernels);

    % Input validation
    if ~iscell(kernels)
        error('Input must be a cell array of kernel filenames.');
    end

    % Iterate and unload each kernel
    for i = 1:numel(kernels)
        kernelName = kernels{i};
        % Locate the file, assuming it was loaded using its full path found by WHICH
        kernelPath = which(kernelName);

        if isempty(kernelPath)
            % If 'which' doesn't find it, it might have been loaded by name
            % or is simply not on the path/not loaded.
            % We'll try to unload by the provided kernelName as a fallback.
            warning('SPICE kernel "%s" not found on MATLAB path. Attempting to unload by name.', kernelName);
            kernelToUnload = kernelName;
        else
            % If 'which' found it, use the full path for unloading.
            kernelToUnload = kernelPath;
        end

        % Attempt to unload
        try
            cspice_unload(kernelToUnload);
        catch ME
            % cspice_unload will error if the kernel wasn't loaded
            % or if there's another issue.
            warning('Failed to unload SPICE kernel "%s": %s', kernelToUnload, ME.message);
        end
    end
end