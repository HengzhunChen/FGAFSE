function plot_fractional_dispersion_1d()
%PLOT_FRACTIONAL_DISPERSION_1D Compare fractional and quadratic dispersion.
%   Plot the FGA probability density for several 1 < alpha < 2 values and
%   for the standard alpha = 2 Schrodinger equation.  

% Add project source paths without depending on the caller's current folder.
thisFile = mfilename('fullpath');
thisDir = fileparts(thisFile);
projectRoot = fileparts(fileparts(thisDir));
addpath(projectRoot);
FGAFSE_startup();

% ------------------------------------------------------------
% Parameters
% ------------------------------------------------------------
alphaList = [1.2, 1.6, 2.0];
vepsExp = -8;
veps = 2^vepsExp;
delta = veps;
rightX = 8;
center = rightX / 2;
finalTime = 5;
solver = 'rk4';
densityFloor = 1e-18;

outputFolder = fullfile(thisDir, 'figures');
if ~exist(outputFolder, 'dir')
    mkdir(outputFolder);
end

initWave = @(x, epsilon) initial_wave(x, epsilon, center);
numAlpha = length(alphaList);
density = cell(numAlpha, 1);

fprintf('One-dimensional fractional-dispersion demo\n');
fprintf('epsilon = 2^(%d) = %.6e, delta = %.6e, T = %.3f\n', ...
    vepsExp, veps, delta, finalTime);

% ------------------------------------------------------------
% FGA solutions
% ------------------------------------------------------------
for k = 1 : numAlpha
    alpha = alphaList(k);
    [wave, x] = FGA1d(alpha, vepsExp, finalTime, rightX, ...
        initWave, @potential, delta, solver);
    density{k} = abs(wave).^2;

    probabilityMass = sum(density{k}) * veps;
    fprintf('alpha = %.1f: integral |psi|^2 dx = %.8f\n', ...
        alpha, probabilityMass);
end

% ------------------------------------------------------------
% Plotting
% ------------------------------------------------------------

% Use solid colorblind-friendly curves for fractional orders and reserve a
% purple dashed curve for the quadratic reference case.
fractionalColors = [0.0000, 0.4470, 0.7410; ...
                    0.8500, 0.3250, 0.0980; ...
                    0.9290, 0.6940, 0.1250; ...
                    0.4940, 0.1840, 0.5560; ...
                    0.4660, 0.6740, 0.1880];

fig = figure('Color', 'w', 'Position', [100, 100, 1500, 430]);

% Linear scale: show the bulk of each probability distribution.
subplot(1, 3, 1)
hold on
for k = 1 : numAlpha
    plot_density_curve(x, density{k}, alphaList(k), k, numAlpha, ...
        fractionalColors);
end
hold off
box on
grid on
xlim([0, rightX])
xlabel('x')
ylabel('|\psi_{FGA}(x)|^2')
title('Probability density')
legend('Location', 'best')

% Semilog scale: reveal small spatial tails on both sides of the packet.
subplot(1, 3, 2)
hold on
for k = 1 : numAlpha
    visibleDensity = density{k};
    visibleDensity(visibleDensity < densityFloor) = NaN;
    plot_density_curve(x, visibleDensity, alphaList(k), k, numAlpha, ...
        fractionalColors);
end
hold off
box on
grid on
set(gca, 'YScale', 'log')
xlim([0, rightX])
ylim([densityFloor, 1e2])
xlabel('x')
ylabel('|\psi_{FGA}(x)|^2')
title('Spatial tails')

% Log-log scale: algebraic tails appear approximately straight, whereas
% the alpha = 2 Gaussian tail bends rapidly downward.
subplot(1, 3, 3)
hold on
tailStart = 4 * sqrt(veps);

for k = 1 : numAlpha
    [~, peakIndex] = max(density{k});
    numOffsets = min(peakIndex - 1, length(x) - peakIndex);
    offsets = (1 : numOffsets)';
    radius = offsets * veps;
    symmetricTail = 0.5 * (density{k}(peakIndex - offsets) ...
        + density{k}(peakIndex + offsets));
    visible = radius >= tailStart & symmetricTail >= densityFloor;
    plot_density_curve(radius(visible), symmetricTail(visible), ...
        alphaList(k), k, numAlpha, fractionalColors);
end
hold off
box on
grid on
set(gca, 'XScale', 'log', 'YScale', 'log')
xlabel('distance from density peak')
ylabel('Symmetrized probability density')
title('Tail decay')

% Typography and axes styling.
axesHandles = findall(fig, 'Type', 'axes');
for k = 1 : numel(axesHandles)
    ax = axesHandles(k);
    set(ax, 'FontSize', 13, 'LineWidth', 1.1, ...
        'TickDir', 'out', 'TickLength', [0.015, 0.015]);
    ax.XLabel.FontSize = 15;
    ax.YLabel.FontSize = 15;
    ax.Title.FontSize = 15;
end
legendHandles = findall(fig, 'Type', 'legend');
set(legendHandles, 'FontSize', 13, 'Box', 'on');

print(fig, fullfile(outputFolder, 'fractional_dispersion_1d.png'), ...
    '-dpng', '-r300');
print(fig, fullfile(outputFolder, 'fractional_dispersion_1d.eps'), ...
    '-depsc');

end


function plot_density_curve(x, y, alpha, index, numAlpha, fractionalColors)
%PLOT_DENSITY_CURVE Apply consistent styling to fractional/reference curves.

    if index == numAlpha
        plot(x, y, '--', 'Color', [0.4940, 0.1840, 0.5560], ...
            'LineWidth', 2.2, ...
            'DisplayName', '\alpha = 2');
    else
        colorIndex = mod(index - 1, size(fractionalColors, 1)) + 1;
        plot(x, y, 'Color', fractionalColors(colorIndex, :), ...
            'LineStyle', '-', 'LineWidth', 1.8, ...
            'DisplayName', sprintf('\\alpha = %.1f', alpha));
    end
end

function u0 = initial_wave(x, veps, center)
    beta = 1;
    u0 = (128 / pi)^(1/4) * exp(-64 * (x - center).^2) ...
        .* exp(1i / veps * beta * x);
end

function [V, DV, D2V] = potential(Q)
    a = 0.1;
    center = 4;
    V = a * (Q - center).^2;
    DV = 2 * a * (Q - center);
    D2V = 2 * a * ones(size(Q));
end
