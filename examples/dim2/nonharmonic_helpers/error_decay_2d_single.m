function result = error_decay_2d_single(alpha, vepsExp, right_x, finalTime, ...
        initWave, potential, useCachedTSSA, deltaPower, fgaDt, tssaFolder)
%ERROR_DECAY_2D_SINGLE Compute one alpha without plotting or shared output.

if nargin < 7 || isempty(useCachedTSSA)
    useCachedTSSA = false;
end
if nargin < 8 || isempty(deltaPower)
    deltaPower = 1;
end
if nargin < 9
    fgaDt = [];
end
if nargin < 10 || isempty(tssaFolder)
    tssaFolder = fullfile(pwd, 'data', 'tssa_2d');
end
if ~isscalar(alpha) || ~isnumeric(alpha) || ~isfinite(alpha)
    error('FGAFSE:InvalidAlpha', 'alpha must be a finite scalar.');
end
if ~isnumeric(vepsExp) || ~isvector(vepsExp)
    error('FGAFSE:InvalidVepsExp', 'vepsExp must be a numeric vector.');
end
if ~exist(tssaFolder, 'dir')
    mkdir(tssaFolder);
end

vepsExp = vepsExp(:).';
nveps = numel(vepsExp);
err_L2 = zeros(1, nveps);
time_FGA = zeros(1, nveps);
time_TSSA = zeros(1, nveps);
solver = 'rk4';

for j = 1:nveps
    veps = 2 ^ vepsExp(j);
    dx = veps;
    fprintf('alpha = %g, veps = %e\n', alpha, veps);

    timer_FGA = tic;
    delta = veps ^ deltaPower;
    [w, xx_w] = FGA2d(alpha, vepsExp(j), finalTime, right_x, ...
        initWave, potential, delta, solver, fgaDt);
    time_FGA(j) = toc(timer_FGA);
    fprintf('FGA2d time: %f seconds\n', time_FGA(j));

    dt = veps ^ 2;
    % dt = veps;
    cacheFile = fullfile(tssaFolder, sprintf( ...
        'alpha_%.12g_vepsExp_%d_T_%0.17g_R_%0.17g.mat', ...
        alpha, vepsExp(j), finalTime, right_x));
    legacyCacheFile = fullfile(tssaFolder, sprintf( ...
        'alpha_%0.17g_vepsExp_%d_T_%0.17g_R_%0.17g.mat', ...
        alpha, vepsExp(j), finalTime, right_x));
    timer_TSSA = tic;
    if useCachedTSSA && exist(cacheFile, 'file')
        cached = load(cacheFile, 'u', 'xx_u');
        u = cached.u;
        xx_u = cached.xx_u;
        fprintf('Loaded TSSA2d result: %s\n', cacheFile);
    elseif useCachedTSSA && exist(legacyCacheFile, 'file')
        cached = load(legacyCacheFile, 'u', 'xx_u');
        u = cached.u;
        xx_u = cached.xx_u;
        fprintf('Loaded legacy TSSA2d result: %s\n', legacyCacheFile);
    else
        [u, xx_u] = TSSA2d(alpha, vepsExp(j), finalTime, right_x, dt, ...
            initWave, potential);
        save(cacheFile, 'u', 'xx_u');
        fprintf('Saved TSSA2d result: %s\n', cacheFile);
    end
    time_TSSA(j) = toc(timer_TSSA);
    fprintf('TSSA2d result time: %f seconds\n', time_TSSA(j));

    if ~isequal(size(u), size(w)) || ~isequal(size(xx_u), size(xx_w))
        error('FGAFSE:MeshMismatch', ...
            'FGA2d and TSSA2d returned incompatible meshes.');
    end
    err_L2(j) = sqrt(sum(sum(abs(u - w).^2)) * dx * dx);
    fprintf('L2 distance between FGA and TSSA: %e\n', err_L2(j));
end

result = struct('alpha', alpha, 'vepsExp', vepsExp, ...
    'err_L2', err_L2, 'time_FGA', time_FGA, 'time_TSSA', time_TSSA, ...
    'right_x', right_x, 'final_time', finalTime, ...
    'deltaPower', deltaPower, 'fgaDt', fgaDt);
end
