%% Wave evolution demo for 2D TSSA and FGA
%
% This script saves separate density GIF animations and PNG snapshots for
% TSSA and FGA. It also reports wave mass and boundary mass, but intentionally
% omits the Q and P trajectory diagnostics used by the 1D demo.

% Add functions into file path
testPath = fileparts(mfilename('fullpath'));
cd ../
FGAFSE_startup();
cd ./test


% ------------------------------------------------------------
% Parameters to edit
% ------------------------------------------------------------
alpha = 1.1;
vepsExp = -6;
right_x = 8;
final_time = 6;
deltaPower = 1;

veps = 2^vepsExp;
dx = veps;
dtTSSA = veps^2;
dtFGA = 1e-2;
delta = veps^deltaPower;
frame_dt = 0.05;
ntTSSA = floor(final_time / dtTSSA + 1e-6);
actualTSSAFinalTime = ntTSSA * dtTSSA;

saveGif = true;
savePng = true;
gifDelayTime = 0.15;
boundaryWidth = min(0.25, right_x / 10);

initWave = @initWavefun;
potential = @potentialfun;

figureFolder = fullfile(testPath, 'figures', 'wave_evolution_2d');
dataFolder = fullfile(testPath, 'data', 'wave_evolution_2d');
if ~exist(figureFolder, 'dir')
    mkdir(figureFolder);
end
if ~exist(dataFolder, 'dir')
    mkdir(dataFolder);
end

fprintf('2D wave evolution demo\n');
fprintf('alpha = %.6f, veps = 2^(%d) = %.6e\n', alpha, vepsExp, veps);
fprintf('right_x = %.6f, final_time = %.6f\n', right_x, final_time);
fprintf('delta = veps^(%.6f) = %.6e\n', deltaPower, delta);
fprintf('actual TSSA final time = %.6f\n', actualTSSAFinalTime);
fprintf('dx = %.6e, TSSA dt = %.6e, FGA dt = %.6e, frame_dt = %.6e\n\n', ...
        dx, dtTSSA, dtFGA, frame_dt);

% ------------------------------------------------------------
% Frame time setup
% ------------------------------------------------------------

requestedTimes = 0 : frame_dt : actualTSSAFinalTime;
if requestedTimes(end) < actualTSSAFinalTime
    requestedTimes = [requestedTimes, actualTSSAFinalTime];
end
[frameSteps, uniqueFrameIdx] = unique(round(requestedTimes / dtTSSA), 'stable');
validFrameIdx = frameSteps >= 0 & frameSteps <= ntTSSA;
frameSteps = frameSteps(validFrameIdx);
uniqueFrameIdx = uniqueFrameIdx(validFrameIdx);
frameTimesTSSA = frameSteps * dtTSSA;
frameTimesFGA = requestedTimes(uniqueFrameIdx);
nFrames = length(frameTimesTSSA);

fprintf('Number of frames: %d\n', nFrames);
fprintf('First frame time = %.6f, last frame time = %.6f\n\n', ...
        frameTimesFGA(1), frameTimesFGA(end));

% ------------------------------------------------------------
% Compute wave evolution frames
% ------------------------------------------------------------

fprintf('Computing TSSA frames...\n');
[uTSSA, xxTSSA] = TSSA_frames(alpha, vepsExp, final_time, right_x, ...
                              dtTSSA, frameSteps, initWave, potential);

fprintf('\nComputing FGA frames...\n');
[uFGA, xxFGA] = FGA_frames(alpha, vepsExp, frameTimesFGA, right_x, ...
                           dtFGA, delta, initWave, potential);

if ~isequal(size(xxTSSA), size(xxFGA)) || ...
        max(abs(xxTSSA(:) - xxFGA(:))) > 10 * eps
    error('TSSA and FGA grids are not aligned.');
end

% ------------------------------------------------------------
% Wave diagnostics
% ------------------------------------------------------------

tssaMass = wave_mass(uTSSA, dx);
fgaMass = wave_mass(uFGA, dx);
tssaBoundaryMass = boundary_mass(uTSSA, xxTSSA, dx, right_x, boundaryWidth);
fgaBoundaryMass = boundary_mass(uFGA, xxFGA, dx, right_x, boundaryWidth);

fprintf('\nTSSA mass drift: %.6e\n', ...
        max(abs(tssaMass - tssaMass(1))) / tssaMass(1));
fprintf('FGA mass drift:  %.6e\n', ...
        max(abs(fgaMass - fgaMass(1))) / fgaMass(1));
fprintf('Max TSSA boundary mass ratio: %.6e\n', ...
        max(tssaBoundaryMass ./ tssaMass));
fprintf('Max FGA boundary mass ratio:  %.6e\n\n', ...
        max(fgaBoundaryMass ./ fgaMass));

% ------------------------------------------------------------
% Save wave GIFs and PNG frames
% ------------------------------------------------------------

densityMax = max([abs(uTSSA(:)).^2; abs(uFGA(:)).^2]);
if densityMax == 0
    densityMax = 1;
end

if saveGif || savePng
    fprintf('Saving TSSA GIF/PNG frames...\n');
    save_wave_outputs(xxTSSA, uTSSA, frameTimesTSSA, alpha, veps, ...
                      tssaMass, tssaBoundaryMass, right_x, densityMax, ...
                      'TSSA', figureFolder, 'tssa', saveGif, savePng, ...
                      gifDelayTime);

    fprintf('Saving FGA GIF/PNG frames...\n');
    save_wave_outputs(xxFGA, uFGA, frameTimesFGA, alpha, veps, ...
                      fgaMass, fgaBoundaryMass, right_x, densityMax, ...
                      'FGA', figureFolder, 'fga', saveGif, savePng, ...
                      gifDelayTime);
end

% ------------------------------------------------------------
% Save numerical data
% ------------------------------------------------------------

save(fullfile(dataFolder, 'wave_evolution_data.mat'), ...
     'alpha', 'vepsExp', 'veps', 'right_x', 'final_time', ...
     'actualTSSAFinalTime', 'deltaPower', 'delta', 'dx', 'dtTSSA', ...
     'dtFGA', 'frame_dt', 'frameTimesTSSA', 'frameTimesFGA', ...
     'xxTSSA', 'uTSSA', 'xxFGA', 'uFGA', 'tssaMass', 'fgaMass', ...
     'tssaBoundaryMass', 'fgaBoundaryMass', 'boundaryWidth');

fprintf('Done. Figure folder:\n%s\n', figureFolder);
fprintf('Data folder:\n%s\n', dataFolder);

% ------------------------------------------------------------
% Auxiliary functions: frame generation
% ------------------------------------------------------------

function [uFrames, xx] = TSSA_frames(alpha, vepsExp, finalTime, right_x, ...
                                     dt, frameSteps, initWavefun, potential)
    veps = 2^vepsExp;
    dx = veps;
    nx = floor(right_x / dx);
    nt = floor(finalTime / dt + 1e-6);
    maxFrameStep = min(max(frameSteps), nt);

    if mod(nx, 2) ~= 0
        error('TSSA_frames expects an even number of Fourier modes. nx = %d.', nx);
    end

    x = linspace(0, right_x, nx + 1);
    x = x(1 : end - 1)';
    xx = x * ones(1, nx);

    [V, ~, ~] = potential(xx, xx');
    u = initWavefun(xx, xx', veps);

    p = [0 : nx/2 - 1, -nx/2 : -1] * 2 * pi / right_x;
    pp = p' * ones(1, nx);
    fLaplace = sqrt(pp.^2 + pp'.^2).^alpha / alpha;
    kineticPhase = exp(-1i * veps^(alpha - 1) * fLaplace * dt);
    potentialHalfPhase = exp(-1i / veps * V * dt / 2);

    uFrames = zeros(nx, nx, length(frameSteps));
    frameIndex = 1;

    if frameSteps(frameIndex) == 0
        uFrames(:, :, frameIndex) = u;
        frameIndex = frameIndex + 1;
    end

    for step = 1 : maxFrameStep
        u = potentialHalfPhase .* u;
        u = fft2(u);
        u = kineticPhase .* u;
        u = ifft2(u);
        u = potentialHalfPhase .* u;

        while frameIndex <= length(frameSteps) && frameSteps(frameIndex) == step
            uFrames(:, :, frameIndex) = u;
            frameIndex = frameIndex + 1;
        end
    end
end

function [uFrames, xx] = FGA_frames(alpha, vepsExp, frameTimes, right_x, ...
                                    dt, delta, initWavefun, potential)
    veps = 2^vepsExp;
    dx = veps;
    nx = floor(right_x / dx);
    dy = dx;
    ny = nx;
    nydq = floor(2^(-vepsExp / 2) / 2);
    kernelSize = floor(2^(-vepsExp / 2)) * 8;

    x = linspace(0, right_x, nx + 1);
    x = x(1 : end - 1)';
    xx = x * ones(1, nx);
    u0 = initWavefun(xx, xx', veps);

    [A, S, Q, P, nGB] = initial_decomposition_2d( ...
        u0, veps, dy, ny, kernelSize, nydq);
    DzQ = repmat([1, 0, 0, 1], nGB, 1);
    DzP = repmat(-1i * [1, 0, 0, 1], nGB, 1);
    odes = @(QValue, PValue, DzQValue, DzPValue) ...
        odes_delta_2d(QValue, PValue, DzQValue, DzPValue, ...
                      alpha, delta, potential);

    uFrames = zeros(nx, nx, length(frameTimes));
    previousTime = 0;

    for frameIndex = 1 : length(frameTimes)
        interval = frameTimes(frameIndex) - previousTime;
        if interval > 0
            [A, S, Q, P, DzQ, DzP] = time_evolution( ...
                A, S, Q, P, DzQ, DzP, dt, interval, odes);
        end

        fprintf('  FGA frame %d/%d, time = %.6f\n', ...
                frameIndex, length(frameTimes), frameTimes(frameIndex));
        uFrames(:, :, frameIndex) = wave_reconstruction_2d( ...
            veps, A, S, Q, P, nGB, dx, nx, kernelSize);
        previousTime = frameTimes(frameIndex);
    end
end

% ------------------------------------------------------------
% Auxiliary functions: diagnostics and output
% ------------------------------------------------------------

function massValues = wave_mass(uFrames, dx)
    massValues = zeros(size(uFrames, 3), 1);
    for frameIndex = 1 : size(uFrames, 3)
        density = abs(uFrames(:, :, frameIndex)).^2;
        massValues(frameIndex) = sum(density(:)) * dx^2;
    end
end

function massValues = boundary_mass(uFrames, xx, dx, right_x, boundaryWidth)
    boundaryMask = xx < boundaryWidth | xx > right_x - boundaryWidth | ...
                   xx' < boundaryWidth | xx' > right_x - boundaryWidth;
    massValues = zeros(size(uFrames, 3), 1);
    for frameIndex = 1 : size(uFrames, 3)
        density = abs(uFrames(:, :, frameIndex)).^2;
        massValues(frameIndex) = sum(density(boundaryMask)) * dx^2;
    end
end

function save_wave_outputs(xx, uFrames, frameTimes, alpha, veps, massValues, ...
                           boundaryMassValues, right_x, densityMax, methodName, ...
                           figureFolder, filePrefix, saveGif, savePng, ...
                           gifDelayTime)
    gifFile = fullfile(figureFolder, [filePrefix, '_wave_evolution.gif']);
    frameFolder = fullfile(figureFolder, 'frames', filePrefix);
    if savePng && ~exist(frameFolder, 'dir')
        mkdir(frameFolder);
    end

    fig = figure('Visible', 'off', 'Color', 'w', ...
                 'Position', [100, 100, 850, 700]);

    for frameIndex = 1 : length(frameTimes)
        clf(fig);
        density = abs(uFrames(:, :, frameIndex)).^2;
        boundaryRatio = boundaryMassValues(frameIndex) / massValues(frameIndex);

        contourf(xx, xx', density, 24, 'LineColor', 'none');
        colormap(fig, 'jet');
        colorbar;
        caxis([0, densityMax]);
        axis equal;
        xlim([0, right_x]);
        ylim([0, right_x]);
        xlabel('x_1');
        ylabel('x_2');
        title(sprintf(['%s density: t = %.6f, alpha = %.3f, eps = %.3e, ' ...
                       'mass = %.6e, boundary/mass = %.3e'], ...
                      methodName, frameTimes(frameIndex), alpha, veps, ...
                      massValues(frameIndex), boundaryRatio));
        drawnow;

        if savePng
            pngFile = fullfile(frameFolder, ...
                sprintf('%s_frame_%04d.png', filePrefix, frameIndex - 1));
        else
            pngFile = [tempname, '.png'];
        end
        saveas(fig, pngFile, 'png');

        if saveGif
            frameImage = imread(pngFile);
            [imind, cm] = rgb2ind(frameImage, 256);
            if frameIndex == 1
                imwrite(imind, cm, gifFile, 'gif', 'Loopcount', inf, ...
                        'DelayTime', gifDelayTime);
            else
                imwrite(imind, cm, gifFile, 'gif', 'WriteMode', 'append', ...
                        'DelayTime', gifDelayTime);
            end
        end

        if ~savePng
            delete(pngFile);
        end
    end

    close(fig);
end

% ------------------------------------------------------------
% Auxiliary functions: test problem
% ------------------------------------------------------------

function u0 = initWavefun(X, Y, veps)
    r = sqrt( (X - 4).^2 + (Y - 4).^2 );
    u0 = sqrt(128 / pi) * exp(-64 * r.^2) ...
        .* exp(1i * (Y - 1) / veps);
end

function [V, DV, D2V] = potentialfun(Q1, Q2)
    a = 0.1;
    b1 = 4; b2 = 4;
    V = a * ((Q1 - b1).^2 + (Q2 - b2).^2);
    DV = 2 * a * [Q1 - b1, Q2 - b2];
    D2V = repmat(2 * a * [1, 0, 0, 1], size(Q1));
end
