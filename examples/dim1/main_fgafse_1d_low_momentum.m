%% 1D FGA low-momentum / weak-force stress test
%
% Compare p0=0.2 and p0=0 at the double-well minimum q0=5.  Full complex
% wavefunction errors are measured against an unregularized Strang reference.

exampleDir = fileparts(mfilename('fullpath'));
cd(fullfile(exampleDir, '..', '..'));
FGAFSE_startup();
cd(exampleDir);

alpha = 1.7;
vepsExp = [-6, -7, -8, -9, -10];
right_x = 8;
final_time = 1;
deltaPower = 1;
fgaDt = 5e-4; % Change this value to study FGA ODE time-step sensitivity.
p0Values = [0.2, 0];
dxLabel = 'dx=epsilon';
dtLabel = 'dt=epsilon^2';

err_L2 = zeros(numel(p0Values), numel(vepsExp));
time_FGA = zeros(size(err_L2));
time_TSSA = zeros(size(err_L2));
for ip = 1:numel(p0Values)
    p0 = p0Values(ip);
    fprintf('\nLow-momentum study: p0 = %g\n', p0);
    for j = 1:numel(vepsExp)
        veps = 2^vepsExp(j);
        dx = veps;
        dt = veps^2;
        fprintf('alpha = %.2f, p0 = %.3g, epsilon = %.6e\n', alpha, p0, veps);

        timer = tic;
        [w, xw] = FGA1d(alpha, vepsExp(j), final_time, right_x, ...
            @(x, e) initWave(x, e, p0), @potential, veps^deltaPower, 'rk4', fgaDt);
        time_FGA(ip, j) = toc(timer);

        timer = tic;
        [u, xu] = TSSA1d(alpha, vepsExp(j), final_time, right_x, dt, dx, ...
            @(x, e) initWave(x, e, p0), @potential);
        time_TSSA(ip, j) = toc(timer);
        if numel(xw) ~= numel(xu) || max(abs(xw - xu)) > 10*eps
            error('FGA and TSSA grids are not aligned.');
        end
        err_L2(ip, j) = sqrt(sum(abs(u - w).^2) * dx);
        fprintf('L2 error = %.6e\n', err_L2(ip, j));
    end
end

% Save all values needed to reproduce the plot and its discretization.
save('L2_error_low_momentum.mat', 'alpha', 'vepsExp', 'right_x', ...
    'final_time', 'deltaPower', 'fgaDt', 'p0Values', 'err_L2', 'time_FGA', ...
    'time_TSSA', 'dxLabel', 'dtLabel');

figure('Name', '1D low-momentum FGA stress test');
hold on; box on; grid on;
styles = {'-*', '-o'};
leg = cell(1, numel(p0Values));
slopes = zeros(size(p0Values));
for ip = 1:numel(p0Values)
    plot(-vepsExp, log2(err_L2(ip, :)), styles{ip}, 'LineWidth', 1);
    fitCoeff = polyfit(vepsExp, log2(err_L2(ip, :)), 1);
    slopes(ip) = fitCoeff(1);
    leg{ip} = sprintf('$p_0=%g$, slope: %.3f', p0Values(ip), slopes(ip));
end
xlabel('$-\log_2(\varepsilon)$', 'Interpreter', 'latex', 'FontSize', 15);
ylabel('$\log_2$($L^2$ error)', 'Interpreter', 'latex', 'FontSize', 15);
legend(leg, 'Interpreter', 'latex', 'Location', 'best');
set(gca, 'FontSize', 14, 'Box', 'on', 'TickDir', 'out');
saveas(gcf, 'L2err_low_momentum.eps', 'epsc');
saveas(gcf, 'L2err_low_momentum.png', 'png');
fprintf('\nFitted rates: p0=0.2: %.3f, p0=0: %.3f\n', slopes(1), slopes(2));
fprintf('Reference refinement: %s, %s; no phase alignment or renormalization.\n', ...
    dxLabel, dtLabel);

function u0 = initWave(x, veps, p0)
    q0 = 5;
    u0 = (128 / pi)^(1/4) * exp(-64 * (x - q0).^2) ...
        .* exp(1i * p0 * (x - q0) / veps);
end

function [V, DV, D2V] = potential(Q)
    y = Q - 4;
    V = (y.^2 - 1).^2 / 16;
    DV = y .* (y.^2 - 1) / 4;
    D2V = (3 * y.^2 - 1) / 4;
end
