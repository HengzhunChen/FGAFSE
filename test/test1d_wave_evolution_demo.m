%% Wave evolution demo for 1D TSSA and FGA
%
% This script saves separate GIF animations and PNG snapshots for the TSSA
% and FGA wave functions. It is intended as a diagnostic for checking whether
% the wave packet reaches the boundary and re-enters because of the periodic
% boundary condition in the Fourier spectral TSSA solver.
%
% Note on frame generation:
% The solvers in src/dim1 return only the final wave at a requested time.
% For TSSA, this script repeats the Strang splitting loop locally so that
% all animation frames can be collected in one time march instead of
% recomputing from t = 0 for every frame. For FGA, this script mirrors the
% FGA1d workflow locally because the animation diagnostics need the hidden
% trajectory variables Q and P, which FGA1d does not return.

% Add functions into file path
testPath = fileparts(mfilename('fullpath'));
cd ../
FGAFSE_startup();
cd ./test


% ------------------------------------------------------------
% Parameters to edit
% ------------------------------------------------------------
alpha = 1.1;
vepsExp = -9;
right_x = 3;
final_time = 4;

veps = 2^vepsExp;
dx = veps;
dt = veps^2;
frame_dt = 0.05;
ntTSSA = floor(final_time / dt + 1e-6);
actualTSSAFinalTime = ntTSSA * dt;

saveGif = true;
savePng = true;
gifDelayTime = 0.15;

boundaryWidth = min(0.25, right_x / 10);
pSmallThreshold = sqrt(veps);

initWave = @initWavefun;
potential = @potentialfun;

figureFolder = fullfile(testPath, 'figures', 'wave_evolution_1d');
dataFolder = fullfile(testPath, 'data', 'wave_evolution_1d');
if ~exist(figureFolder, 'dir')
    mkdir(figureFolder);
end
if ~exist(dataFolder, 'dir')
    mkdir(dataFolder);
end

fprintf('1D wave evolution demo\n');
fprintf('alpha = %.6f, veps = 2^(%d) = %.6e\n', alpha, vepsExp, veps);
fprintf('right_x = %.6f, final_time = %.6f\n', right_x, final_time);
fprintf('actual TSSA final time = %.6f\n', actualTSSAFinalTime);
fprintf('dx = %.6e, dt = %.6e, frame_dt = %.6e\n\n', dx, dt, frame_dt);

% ------------------------------------------------------------
% Frame time setup
% ------------------------------------------------------------

% The requested frame times are rounded to TSSA time steps. This keeps the
% displayed TSSA states exactly on the computed time grid.
% NOTE: frame_dt is different with dt for TSSA and FGA.
requestedTimes = 0 : frame_dt : actualTSSAFinalTime;
if requestedTimes(end) < actualTSSAFinalTime
    requestedTimes = [requestedTimes, actualTSSAFinalTime];
end
[frameSteps, uniqueFrameIdx] = unique(round(requestedTimes / dt), 'stable');
validFrameIdx = frameSteps >= 0 & frameSteps <= ntTSSA;
frameSteps = frameSteps(validFrameIdx);
uniqueFrameIdx = uniqueFrameIdx(validFrameIdx);
frameTimesTSSA = frameSteps * dt;
frameTimesFGA = requestedTimes(uniqueFrameIdx);
nFrames = length(frameTimesTSSA);

fprintf('Number of frames: %d\n', nFrames);
fprintf('First TSSA frame time = %.6f, last TSSA frame time = %.6f\n', ...
        frameTimesTSSA(1), frameTimesTSSA(end));
fprintf('First FGA frame time = %.6f, last FGA frame time = %.6f\n\n', ...
        frameTimesFGA(1), frameTimesFGA(end));

% ------------------------------------------------------------
% Compute wave evolution frames
% ------------------------------------------------------------

fprintf('Computing TSSA frames...\n');
[uTSSA, xTSSA] = TSSA_frames(alpha, vepsExp, final_time, right_x, dt, dx, ...
                            frameSteps, initWave, potential);

fprintf('\nComputing FGA frames...\n');
[uFGA, xFGA, fgaQPStats] = FGA_frames(alpha, vepsExp, frameTimesFGA, ...
                                     right_x, pSmallThreshold, initWave, ...
                                     potential);

if length(xTSSA) ~= length(xFGA) || max(abs(xTSSA - xFGA)) > 10 * eps
    error('TSSA and FGA grids are not aligned.');
end

% ------------------------------------------------------------
% Mass, boundary, and trajectory diagnostics
% ------------------------------------------------------------

tssaMass = wave_mass(uTSSA, dx);
fgaMass = wave_mass(uFGA, dx);
tssaBoundaryMass = boundary_mass(uTSSA, xTSSA, dx, right_x, boundaryWidth);
fgaBoundaryMass = boundary_mass(uFGA, xFGA, dx, right_x, boundaryWidth);
fgaMinQ = [fgaQPStats.minQ]';
fgaMaxQ = [fgaQPStats.maxQ]';
fgaNumQOutside = [fgaQPStats.numQOutside]';
fgaNumQLeft = [fgaQPStats.numQLeft]';
fgaNumQRight = [fgaQPStats.numQRight]';
fgaMinP = [fgaQPStats.minP]';
fgaMaxP = [fgaQPStats.maxP]';
fgaMinAbsP = [fgaQPStats.minAbsP]';
fgaNumPSmall = [fgaQPStats.numPSmall]';
fgaNumGB = [fgaQPStats.nGB]';

fprintf('\nTSSA mass drift: %.6e\n', max(abs(tssaMass - tssaMass(1))) / tssaMass(1));
fprintf('FGA mass drift:  %.6e\n', max(abs(fgaMass - fgaMass(1))) / fgaMass(1));
fprintf('Max TSSA boundary mass ratio: %.6e\n', max(tssaBoundaryMass ./ tssaMass));
fprintf('Max FGA boundary mass ratio:  %.6e\n\n', max(fgaBoundaryMass ./ fgaMass));

fprintf('FGA trajectory Q diagnostics\n');
fprintf('  time        min(Q)        max(Q)        #outside    #left    #right    #GB\n');
for frameIndex = 1 : nFrames
    fprintf('  %.6f   %.6e   %.6e   %8d   %6d   %7d   %6d\n', ...
            frameTimesFGA(frameIndex), fgaMinQ(frameIndex), ...
            fgaMaxQ(frameIndex), fgaNumQOutside(frameIndex), ...
            fgaNumQLeft(frameIndex), fgaNumQRight(frameIndex), ...
            fgaNumGB(frameIndex));
end
[maxOutside, maxOutsideIdx] = max(fgaNumQOutside);
fprintf('Max FGA outside Q count: %d / %d at time %.6f\n', ...
        maxOutside, fgaNumGB(maxOutsideIdx), frameTimesFGA(maxOutsideIdx));
fprintf('Overall FGA Q range: [%.6e, %.6e]\n\n', min(fgaMinQ), max(fgaMaxQ));

fprintf('FGA momentum P diagnostics, small threshold = %.6e\n', pSmallThreshold);
fprintf('  time        min(P)        max(P)      min(abs(P))   #small |P|    #GB\n');
for frameIndex = 1 : nFrames
    fprintf('  %.6f   %.6e   %.6e   %.6e   %10d   %6d\n', ...
            frameTimesFGA(frameIndex), fgaMinP(frameIndex), ...
            fgaMaxP(frameIndex), fgaMinAbsP(frameIndex), ...
            fgaNumPSmall(frameIndex), fgaNumGB(frameIndex));
end
[maxSmallP, maxSmallPIdx] = max(fgaNumPSmall);
[minAbsPValue, minAbsPIdx] = min(fgaMinAbsP);
fprintf('Max small-|P| count: %d / %d at time %.6f\n', ...
        maxSmallP, fgaNumGB(maxSmallPIdx), frameTimesFGA(maxSmallPIdx));
fprintf('Minimum abs(P): %.6e at time %.6f\n\n', ...
        minAbsPValue, frameTimesFGA(minAbsPIdx));

% ------------------------------------------------------------
% Save diagnostic figures
% ------------------------------------------------------------

save_fga_qp_diagnostics_plot(frameTimesFGA, fgaMinQ, fgaMaxQ, ...
                             fgaNumQOutside, fgaMinP, fgaMaxP, ...
                             fgaMinAbsP, fgaNumPSmall, right_x, ...
                             pSmallThreshold, figureFolder);

% ------------------------------------------------------------
% Save wave GIFs and PNG frames
% ------------------------------------------------------------

tssaStyle = plot_style(uTSSA);
fgaStyle = plot_style(uFGA);
tssaFrameText = repmat({''}, nFrames, 1);
fgaFrameText = repmat({''}, nFrames, 1);

if saveGif || savePng
    fprintf('Saving TSSA GIF/PNG frames...\n');
    save_wave_outputs(xTSSA, uTSSA, frameTimesTSSA, alpha, veps, ...
                      tssaMass, tssaBoundaryMass, boundaryWidth, right_x, ...
                      tssaStyle, tssaFrameText, 'TSSA', figureFolder, 'tssa', ...
                      saveGif, savePng, gifDelayTime);

    fprintf('Saving FGA GIF/PNG frames...\n');
    save_wave_outputs(xFGA, uFGA, frameTimesFGA, alpha, veps, ...
                      fgaMass, fgaBoundaryMass, boundaryWidth, right_x, ...
                      fgaStyle, fgaFrameText, 'FGA', figureFolder, 'fga', ...
                      saveGif, savePng, gifDelayTime);
end

% ------------------------------------------------------------
% Save numerical data
% ------------------------------------------------------------

save(fullfile(dataFolder, 'wave_evolution_data.mat'), ...
     'alpha', 'vepsExp', 'veps', 'right_x', 'final_time', ...
     'actualTSSAFinalTime', 'dx', 'dt', 'frame_dt', ...
     'frameTimesTSSA', 'frameTimesFGA', 'xTSSA', 'uTSSA', ...
     'xFGA', 'uFGA', 'tssaMass', 'fgaMass', ...
     'tssaBoundaryMass', 'fgaBoundaryMass', 'boundaryWidth', ...
     'fgaQPStats', 'fgaMinQ', 'fgaMaxQ', 'fgaNumQOutside', ...
     'fgaNumQLeft', 'fgaNumQRight', 'fgaMinP', 'fgaMaxP', ...
     'fgaMinAbsP', 'fgaNumPSmall', 'pSmallThreshold', 'fgaNumGB');

fprintf('Done. Figure folder:\n%s\n', figureFolder);
fprintf('Data folder:\n%s\n', dataFolder);

% ------------------------------------------------------------
% Auxiliary functions: frame generation
% ------------------------------------------------------------

function [uFrames, x] = TSSA_frames(alpha, vepsExp, finalTime, right_x, dt, dx, ...
                                    frameSteps, initWavefun, potential)
    veps = 2 ^ vepsExp;

    nx = floor(right_x / dx);
    nt = floor(finalTime / dt + 1e-6);
    maxFrameStep = min(max(frameSteps), nt);

    if mod(nx, 2) ~= 0
        error('TSSA_frames expects an even number of Fourier modes. nx = %d.', nx);
    end

    x = 0 : dx : right_x;
    x = x(1 : end-1)';
    u = initWavefun(x, veps);
    V = potential(x);

    k = [0 : nx/2 - 1, -nx/2 : -1]' * 2 * pi / right_x;
    fLaplace = (abs(k) .^ alpha) / alpha;
    kineticPhase = exp(-1i * (veps ^ (alpha - 1)) * fLaplace * dt);
    potentialHalfPhase = exp(-1i / veps * V * dt / 2);

    uFrames = zeros(nx, length(frameSteps));
    frameIndex = 1;

    if frameSteps(frameIndex) == 0
        uFrames(:, frameIndex) = u;
        frameIndex = frameIndex + 1;
    end

    for step = 1 : maxFrameStep
        u = potentialHalfPhase .* u;
        u = fft(u);
        u = kineticPhase .* u;
        u = ifft(u);
        u = potentialHalfPhase .* u;

        while frameIndex <= length(frameSteps) && frameSteps(frameIndex) == step
            uFrames(:, frameIndex) = u;
            frameIndex = frameIndex + 1;
        end
    end
end

function [uFrames, x, qpStats] = FGA_frames(alpha, vepsExp, frameTimes, right_x, ...
                                           pSmallThreshold, initWavefun, potential)
    uFrames = [];
    x = [];
    qpStats = struct('minQ', {}, 'maxQ', {}, 'numQOutside', {}, ...
                    'numQLeft', {}, 'numQRight', {}, 'minP', {}, ...
                    'maxP', {}, 'minAbsP', {}, 'numPSmall', {}, 'nGB', {});

    for frameIndex = 1 : length(frameTimes)
        fprintf('  FGA frame %d/%d, time = %.6f\n', ...
                frameIndex, length(frameTimes), frameTimes(frameIndex));

        [u, xFrame, qpStats(frameIndex)] = FGA1d_with_QP_stats(alpha, vepsExp, ...
            frameTimes(frameIndex), right_x, pSmallThreshold, initWavefun, ...
            potential);

        if frameIndex == 1
            x = xFrame;
            uFrames = zeros(length(u), length(frameTimes));
        elseif length(xFrame) ~= length(x) || max(abs(xFrame - x)) > 10 * eps
            error('FGA grid changed at frame %d.', frameIndex);
        end

        uFrames(:, frameIndex) = u;
    end
end

function [w, x, qpStats] = FGA1d_with_QP_stats(alpha, vepsExp, finalTime, ...
                                             right_x, pSmallThreshold, ...
                                             initWavefun, potential)
    veps = 2 ^ vepsExp;

    dx = veps;
    nx = floor(right_x / dx);
    dy = dx;
    ny = nx;

    nydq = floor(2 ^ (-vepsExp / 2) / 2);
    kernelSize = floor(2 ^ (-vepsExp / 2)) * 2^5;

    % NOTE: if dt is too large, FGA solution may lose some mass after
    % reflection with potential boundary
    % dt = 1e-2;
    dt = 1e-3;

    x = 0 : dx : right_x;
    x = x(1 : end-1)';
    u0 = initWavefun(x, veps);

    delta = veps;
    odes = @(Q, P, DzQ, DzP) ...
        odes_delta_1d(Q, P, DzQ, DzP, alpha, delta, potential);

    [A0, S0, Q0, P0, nGB] = initial_decomposition_1d(u0, veps, dy, ny, ...
                                                     kernelSize, nydq);
    DzQ0 = ones(size(Q0));
    DzP0 = -1i * ones(size(P0));
    [A, S, Q, P] = time_evolution(A0, S0, Q0, P0, DzQ0, DzP0, dt, finalTime, odes);

    qpStats.minQ = min(Q);
    qpStats.maxQ = max(Q);
    qpStats.numQLeft = sum(Q < 0);
    qpStats.numQRight = sum(Q > right_x);
    qpStats.numQOutside = qpStats.numQLeft + qpStats.numQRight;
    qpStats.minP = min(P);
    qpStats.maxP = max(P);
    qpStats.minAbsP = min(abs(P));
    qpStats.numPSmall = sum(abs(P) < pSmallThreshold);
    qpStats.nGB = nGB;

    w = wave_reconstruction_1d(veps, A, S, Q, P, nGB, x, dx, nx, kernelSize);
end

% ------------------------------------------------------------
% Auxiliary functions: diagnostic plots
% ------------------------------------------------------------

function save_fga_qp_diagnostics_plot(frameTimes, minQ, maxQ, numQOutside, ...
                                      minP, maxP, minAbsP, numPSmall, ...
                                      right_x, pSmallThreshold, outputFolder)
    fig = figure('Visible', 'off', 'Color', 'w', 'Position', [100, 100, 900, 850]);

    subplot(4, 1, 1);
    hold on
    box on
    grid on
    plot(frameTimes, minQ, '-o', 'LineWidth', 1.2);
    plot(frameTimes, maxQ, '-s', 'LineWidth', 1.2);
    plot([frameTimes(1), frameTimes(end)], [0, 0], 'k--');
    plot([frameTimes(1), frameTimes(end)], [right_x, right_x], 'k--');
    hold off
    xlabel('time');
    ylabel('Q range');
    legend('min(Q)', 'max(Q)', 'left boundary', 'right boundary', ...
           'Location', 'best');
    title('FGA trajectory centers');

    subplot(4, 1, 2);
    box on
    grid on
    plot(frameTimes, numQOutside, '-o', 'LineWidth', 1.2);
    xlabel('time');
    ylabel('# outside [0, right\_x]');
    title('FGA trajectory centers outside domain');

    subplot(4, 1, 3);
    hold on
    box on
    grid on
    plot(frameTimes, minP, '-o', 'LineWidth', 1.2);
    plot(frameTimes, maxP, '-s', 'LineWidth', 1.2);
    plot([frameTimes(1), frameTimes(end)], [0, 0], 'k--');
    hold off
    xlabel('time');
    ylabel('P range');
    legend('min(P)', 'max(P)', 'P = 0', 'Location', 'best');
    title('FGA momenta');

    subplot(4, 1, 4);
    hold on
    box on
    grid on
    plot(frameTimes, minAbsP, '-o', 'LineWidth', 1.2);
    plot([frameTimes(1), frameTimes(end)], ...
         [pSmallThreshold, pSmallThreshold], 'k--');
    yyaxis right
    plot(frameTimes, numPSmall, '-s', 'LineWidth', 1.2);
    ylabel('# small |P|');
    yyaxis left
    hold off
    xlabel('time');
    ylabel('min |P|');
    title('Small momentum diagnostic');

    saveas(fig, fullfile(outputFolder, 'fga_qp_diagnostics.png'), 'png');
    close(fig);
end

% ------------------------------------------------------------
% Auxiliary functions: GIF and PNG output
% ------------------------------------------------------------

function save_wave_outputs(x, uFrames, frameTimes, alpha, veps, massValues, ...
                           boundaryMassValues, boundaryWidth, right_x, ...
                           style, frameText, methodName, outputFolder, ...
                           filePrefix, saveGif, savePng, gifDelayTime)
    gifFile = fullfile(outputFolder, [filePrefix, '_wave_evolution.gif']);
    frameFolder = fullfile(outputFolder, 'frames', filePrefix);
    if savePng && ~exist(frameFolder, 'dir')
        mkdir(frameFolder);
    end

    fig = figure('Visible', 'off', 'Color', 'w', 'Position', [100, 100, 900, 650]);

    for frameIndex = 1 : length(frameTimes)
        clf(fig);
        plot_wave_frame(fig, x, uFrames(:, frameIndex), frameTimes(frameIndex), ...
                        alpha, veps, massValues(frameIndex), ...
                        boundaryMassValues(frameIndex), boundaryWidth, ...
                        right_x, style, frameText{frameIndex}, methodName);

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
% Auxiliary functions: plotting utilities
% ------------------------------------------------------------

function plot_wave_frame(fig, x, u, timeValue, alpha, veps, massValue, ...
                         boundaryMassValue, boundaryWidth, right_x, style, ...
                         frameText, methodName)
    density = abs(u) .^ 2;
    realPart = real(u);
    boundaryRatio = boundaryMassValue / massValue;

    ax1 = subplot(2, 1, 1, 'Parent', fig);
    hold(ax1, 'on');
    mark_boundary(ax1, boundaryWidth, right_x, style.densityYLim);
    plot(ax1, x, density, 'b-', 'LineWidth', 1.3);
    hold(ax1, 'off');
    box(ax1, 'on');
    grid(ax1, 'on');
    xlim(ax1, [0, right_x]);
    ylim(ax1, style.densityYLim);
    ylabel(ax1, '|u|^2');
    titleText = sprintf('%s: t = %.6f, alpha = %.3f, eps = %.3e, mass = %.6e, boundary/mass = %.3e', ...
                        methodName, timeValue, alpha, veps, massValue, boundaryRatio);
    if ~isempty(frameText)
        titleText = [titleText, ', ', frameText];
    end
    title(ax1, titleText);

    ax2 = subplot(2, 1, 2, 'Parent', fig);
    hold(ax2, 'on');
    mark_boundary(ax2, boundaryWidth, right_x, style.realYLim);
    plot(ax2, x, realPart, 'r-', 'LineWidth', 1.0);
    hold(ax2, 'off');
    box(ax2, 'on');
    grid(ax2, 'on');
    xlim(ax2, [0, right_x]);
    ylim(ax2, style.realYLim);
    xlabel(ax2, 'x');
    ylabel(ax2, 'Re(u)');
    title(ax2, 'Real part');
end

function mark_boundary(ax, boundaryWidth, right_x, yLim)
    patch(ax, [0, boundaryWidth, boundaryWidth, 0], ...
          [yLim(1), yLim(1), yLim(2), yLim(2)], ...
          [0.92, 0.92, 0.92], 'EdgeColor', 'none', 'FaceAlpha', 0.7);
    patch(ax, [right_x - boundaryWidth, right_x, right_x, right_x - boundaryWidth], ...
          [yLim(1), yLim(1), yLim(2), yLim(2)], ...
          [0.92, 0.92, 0.92], 'EdgeColor', 'none', 'FaceAlpha', 0.7);
end

function style = plot_style(uFrames)
    densityMax = max(abs(uFrames(:)) .^ 2);
    realMax = max(abs(real(uFrames(:))));

    if densityMax == 0
        densityMax = 1;
    end
    if realMax == 0
        realMax = 1;
    end

    style.densityYLim = [0, 1.05 * densityMax];
    style.realYLim = [-1.05 * realMax, 1.05 * realMax];
end

function massValues = wave_mass(uFrames, dx)
    massValues = sum(abs(uFrames) .^ 2, 1)' * dx;
end

function massValues = boundary_mass(uFrames, x, dx, right_x, boundaryWidth)
    idx = x < boundaryWidth | x > right_x - boundaryWidth;
    massValues = sum(abs(uFrames(idx, :)) .^ 2, 1)' * dx;
end

% ------------------------------------------------------------
% Auxiliary functions: test problem
% ------------------------------------------------------------

function u0 = initWavefun(x, veps)
    beta = 1;
    u0 = (128 / pi)^(1/4) * exp(-64 * (x - 1.5).^2) ...
     .* exp(1i / veps * beta * x);
end

function [V, DV, D2V] = potentialfun(Q)
    a = 1;
    b = 1.5;
    V = a * (Q - b).^2;
    DV = 2 * a * (Q - b);
    D2V = 2 * a;
end
