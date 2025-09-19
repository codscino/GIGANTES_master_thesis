function loadSpiceKernels(kernels, isParallel)
%   Load SPICE kernels, optionally in parallel on all workers of a parallel pool.
%
%   Inputs:
%       kernels    - A cell array of kernel filenames (e.g., {'sat441.bsp', 'naif0012.tls'}).
%                    Each kernel file must be discoverable on the MATLAB path
%                    of the session (or worker) attempting to load it.
%       isParallel - (Optional) Logical flag. If true, the function attempts
%                    to load kernels on all workers in a parallel pool.
%                    If false or not provided, kernels are loaded sequentially
%                    in the current MATLAB session. Defaults to false.
%
%   Example:
%       kernels = {'sat441.bsp', 'naif0012.tls', 'pck00010.tpc'};
%       loadSpiceKernels(kernels);          % Loads sequentially
%       loadSpiceKernels(kernels, true);    % Attempts parallel load

    % --- Input validation for 'kernels' ---
    if ~iscell(kernels)
        error('Input ''kernels'' must be a cell array of kernel filenames.');
    end

    % --- Set default for isParallel if not provided ---
    if nargin < 2
        isParallel = false;
    else
        if ~islogical(isParallel) || ~isscalar(isParallel)
            error('Input ''isParallel'' must be a logical scalar.');
        end
    end

    % --- Check for Parallel Computing Toolbox if parallel loading is requested ---
    if isParallel
        if license('test', 'Distrib_Computing_Toolbox') && ~isempty(ver('parallel'))
            pool = gcp('nocreate'); % Get current pool without creating one
            if isempty(pool)
                try
                    parpool(); % Start a new parallel pool
                catch ME
                    warning('Failed to start parallel pool: %s. Proceeding with sequential loading.', ME.message);
                    isParallel = false; % Fallback to sequential if parpool fails
                end
            end
        else
            warning('Parallel Computing Toolbox not available or not installed. Proceeding with sequential loading.');
            isParallel = false; % Fallback to sequential
        end
    end

    % --- Execute kernel loading based on isParallel flag ---
    if isParallel
        % --- Parallel Loading Logic (inside an SPMD block) ---
        spmd
            % This block executes on each worker in the parallel pool
            loadKernelsInternal(kernels, labindex); % Call helper for actual loading
        end
    else
        % --- Sequential Loading Logic (in the current MATLAB session) ---
        loadKernelsInternal(kernels); % Call helper for actual loading
    end
end

% --- Helper function for the actual kernel loading logic ---
% This internal function centralizes the core logic for finding and furnishing kernels.
function loadKernelsInternal(kernels, workerIndex)
    if nargin < 2
        workerIndex = 0; % Use 0 to indicate main thread/sequential
    end

    for i = 1:numel(kernels)
        kernelName = kernels{i};
        % Locate the file on the current MATLAB path (or worker's path)
        kernelPath = which(kernelName);

        if isempty(kernelPath)
            if workerIndex > 0
                warning('Worker %d: SPICE kernel "%s" not found on MATLAB path.', workerIndex, kernelName);
            else
                warning('SPICE kernel "%s" not found on MATLAB path.', kernelName);
            end
            continue;
        end

        % Attempt to load the kernel
        try
            cspice_furnsh(kernelPath);
        catch ME
            if workerIndex > 0
                warning('Worker %d: Failed to furnish SPICE kernel "%s": %s', workerIndex, kernelPath, ME.message);
            else
                warning('Failed to furnish SPICE kernel "%s": %s', kernelPath, ME.message);
            end
        end
    end
end