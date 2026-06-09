%% Space refinement test for the 1D TSSA reference solution
%
% This script checks whether the spatial grid dx = eps is fine enough for
% the time-splitting spectral approximation used as the FGA reference.
%
% The grids are nested:
%
%     dx = eps, eps/2, eps/4, eps/8, ...
%
% A coarse-grid solution is compared with the next finer solution restricted
% to the same grid points. This avoids adding interpolation error to the
% diagnostic.

% Add functions into file path
testPath = fileparts(mfilename('fullpath'));
cd ../
FGAFSE_startup();
cd ./test

% ------------------------------------------------------------
% Parameters to edit
% ------------------------------------------------------------
right_x = 3;
final_time = 2;

alpha_values = 1.1;
vepsExp_values = [-6, -7, -8, -9, -10];

refinementFactors = [1, 2, 4, 8];
timeStepFactor = 1.0;  % For verifying the convergence with dt when necessary

initWave = @initWavefun;
potential = @potentialfun;

folder = fullfile(testPath, 'figures');
if ~exist(folder, 'dir')
    mkdir(folder);
end

fprintf('1D TSSA space refinement test\n');
fprintf('right_x = %.6f, final_time = %.6f\n', right_x, final_time);
fprintf('refinement factors: ');
fprintf('%d ', refinementFactors);
fprintf('\n\n');

% ------------------------------------------------------------
% Main loop
% ------------------------------------------------------------

for ia = 1 : length(alpha_values)
    alpha = alpha_values(ia);

    for ie = 1 : length(vepsExp_values)
        vepsExp = vepsExp_values(ie);
        veps = 2^vepsExp;
        dt = timeStepFactor * veps^2;

        fprintf('alpha = %.6f, veps = 2^(%d) = %.6e, dt = %.6e\n', ...
                alpha, vepsExp, veps, dt);

        dx_values = veps ./ refinementFactors;
        u_values = cell(length(dx_values), 1);
        x_values = cell(length(dx_values), 1);
        mass0 = zeros(length(dx_values), 1);
        massT = zeros(length(dx_values), 1);
        boundaryMax = zeros(length(dx_values), 1);

        for k = 1 : length(dx_values)
            dx = dx_values(k);
            nx = floor(right_x / dx);

            if abs(nx * dx - right_x) > 100 * eps(right_x)
                error('right_x/dx is not an integer for dx = %.16e.', dx);
            end

            if mod(nx, 2) ~= 0
                error('TSSA1d expects an even number of Fourier modes. nx = %d.', nx);
            end

            fprintf('  solving with dx = eps/%d = %.6e, nx = %d\n', ...
                    refinementFactors(k), dx, nx);

            [u_values{k}, x_values{k}] = TSSA1d(alpha, vepsExp, ...
                final_time, right_x, dt, dx, initWave, potential);

            u0 = initWave(x_values{k}, veps);
            mass0(k) = sum(abs(u0) .^ 2) * dx;
            massT(k) = sum(abs(u_values{k}) .^ 2) * dx;
            boundaryMax(k) = boundary_max(u_values{k}, x_values{k}, right_x);
        end

        errToFinest = zeros(length(dx_values), 1);
        relErrToFinest = zeros(length(dx_values), 1);
        errAdjacent = nan(length(dx_values), 1);
        relErrAdjacent = nan(length(dx_values), 1);

        uFine = u_values{end};
        dxFine = dx_values(end);

        for k = 1 : length(dx_values)
            ratio = refinementFactors(end) / refinementFactors(k);
            idx = 1 : ratio : length(uFine);
            uFineOnGrid = uFine(idx);

            errToFinest(k) = sqrt(sum(abs(u_values{k} - uFineOnGrid) .^ 2) * dx_values(k));
            relErrToFinest(k) = errToFinest(k) / ...
                sqrt(sum(abs(uFineOnGrid) .^ 2) * dx_values(k));

            if k < length(dx_values)
                ratioAdjacent = refinementFactors(k + 1) / refinementFactors(k);
                uNextOnGrid = u_values{k + 1}(1 : ratioAdjacent : end);

                errAdjacent(k) = sqrt(sum(abs(u_values{k} - uNextOnGrid) .^ 2) * dx_values(k));
                relErrAdjacent(k) = errAdjacent(k) / ...
                    sqrt(sum(abs(uNextOnGrid) .^ 2) * dx_values(k));
            end
        end

        massDrift = abs(massT - mass0) ./ mass0;

        fprintf('\n');
        fprintf('  factor        dx          L2-to-finest   rel-to-finest   L2-to-next     rel-to-next     mass-drift     boundary-max\n');
        for k = 1 : length(dx_values)
            fprintf('  %6d   %.6e   %.6e   %.6e   %.6e   %.6e   %.6e   %.6e\n', ...
                    refinementFactors(k), dx_values(k), errToFinest(k), ...
                    relErrToFinest(k), errAdjacent(k), relErrAdjacent(k), ...
                    massDrift(k), boundaryMax(k));
        end
        fprintf('\n');

        figure;
        hold on
        box on
        grid on
        loglog(dx_values(1 : end-1), errToFinest(1 : end-1), '-o', ...
               'LineWidth', 1.2);
        loglog(dx_values(1 : end-1), errAdjacent(1 : end-1), '-s', ...
               'LineWidth', 1.2);
        set(gca, 'XDir', 'reverse');
        hold off
        xlabel('dx');
        ylabel('L^2 difference');
        legend('difference to finest grid', 'difference to next finer grid', ...
               'Location', 'best');
        title(['TSSA space refinement, \alpha = ', num2str(alpha), ...
               ', \epsilon = 2^{', num2str(vepsExp), '}']);

        filename = fullfile(folder, ...
            ['test1d_tssa_space_refinement_alpha_', ...
             strrep(num2str(alpha), '.', 'p'), ...
             '_veps_', num2str(-vepsExp), '.png']);
        saveas(gcf, filename, 'png');
    end
end

% ------------------------------------------------------------
% Auxiliary functions
% ------------------------------------------------------------

function u0 = initWavefun(x, veps)
    beta = 1;
    u0 = exp(-64 * (x - 0.5) .^ 2) / sqrt(pi / 64) .* ...
         exp(1i / veps * beta * x);
end

function [V, DV, D2V] = potentialfun(Q)
    % Bad example
    V = 1 + cos(pi * Q);
    DV = -pi * sin(pi * Q);
    D2V = -(pi)^2 * cos(pi * Q);
end

function val = boundary_max(u, x, right_x)
    boundaryWidth = min(0.25, right_x / 10);
    idx = x < boundaryWidth | x > right_x - boundaryWidth;
    val = max(abs(u(idx)));
end
