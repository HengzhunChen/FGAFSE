%% 2D FGA low-momentum / weak-force stress test

exampleDir = fileparts(mfilename('fullpath'));
cd(fullfile(exampleDir, '..', '..'));
FGAFSE_startup();
cd(exampleDir);

alpha = 1.7;
vepsExp = [-5, -6, -7, -8];
right_x = 8;
final_time = 1;
deltaPower = 1;
fgaDt = 1e-3; % FGA ODE step; reduce this to test time-step sensitivity.
quadraticScale = 1/16;
quarticScale = 1/128;
potentialName = sprintf(['V=%.16g*|y|^2 + %.16g*|y|^4, ' ...
    'y=(x1-4,x2-4)'], quadraticScale, quarticScale);
pcValues = [0.2, 0]; % The 2D carrier momentum is (0, pc).

err_L2 = zeros(numel(pcValues), numel(vepsExp));
time_FGA = zeros(size(err_L2));
time_TSSA = zeros(size(err_L2));
for ip = 1:numel(pcValues)
    pc = pcValues(ip);
    for j = 1:numel(vepsExp)
        veps = 2^vepsExp(j);
        dx = veps;
        dt = veps^2;
        fprintf('alpha = %.2f, pc = %g, epsilon = %.6e\n', ...
            alpha, pc, veps);
        timer = tic;
        [w, xw] = FGA2d(alpha, vepsExp(j), final_time, right_x, ...
            @(X,Y,e) initWave(X,Y,e,pc), ...
            @(Q1,Q2) potential(Q1,Q2,quadraticScale,quarticScale), ...
            veps^deltaPower, ...
            'rk4', fgaDt);
        time_FGA(ip,j) = toc(timer);
        timer = tic;
        [u, xu] = TSSA2d(alpha, vepsExp(j), final_time, right_x, dt, ...
            @(X,Y,e) initWave(X,Y,e,pc), ...
            @(Q1,Q2) potential(Q1,Q2,quadraticScale,quarticScale));
        time_TSSA(ip,j) = toc(timer);
        if ~isequal(size(xw), size(xu)) || max(abs(xw(:)-xu(:))) > 10*eps
            error('FGA and TSSA grids are not aligned.');
        end
        err_L2(ip,j) = sqrt(sum(abs(u(:)-w(:)).^2) * dx^2);
        fprintf('L2 error = %.6e\n', err_L2(ip,j));
    end
end

save('L2_error_low_momentum_2d.mat', 'alpha', 'vepsExp', 'right_x', ...
    'final_time', 'deltaPower', 'fgaDt', 'quadraticScale', ...
    'quarticScale', 'potentialName', 'pcValues', 'err_L2', ...
    'time_FGA', 'time_TSSA');

figure('Name', '2D low-momentum FGA stress test');
hold on; box on; grid on;
styles = {'-*', '-o'};
leg = cell(1, numel(pcValues));
slopes = zeros(1, numel(pcValues));
for ip = 1:numel(pcValues)
    plot(-vepsExp, log2(err_L2(ip,:)), styles{ip}, 'LineWidth', 1);
    fitCoeff = polyfit(vepsExp, log2(err_L2(ip,:)), 1);
    slopes(ip) = fitCoeff(1);
    leg{ip} = sprintf('$p_{\\mathrm{c}}=%g$, slope: %.3f', ...
        pcValues(ip), slopes(ip));
end
xlabel('$-\log_2(\varepsilon)$', 'Interpreter', 'latex', 'FontSize', 15);
ylabel('$\log_2$($L^2$ error)', 'Interpreter', 'latex', 'FontSize', 15);
legend(leg, 'Interpreter', 'latex', 'Location', 'best');
set(gca, 'FontSize', 14, 'Box', 'on', 'TickDir', 'out');
saveas(gcf, 'L2err_low_momentum_2d.eps', 'epsc');
saveas(gcf, 'L2err_low_momentum_2d.png', 'png');
fprintf('Fitted rates: moving %.3f, resting %.3f\n', slopes(1), slopes(2));

function u0 = initWave(X, Y, veps, pc)
    q0 = [4, 4];
    u0 = sqrt(128 / pi) * exp(-64 * ((X-q0(1)).^2 + (Y-q0(2)).^2)) ...
        .* exp(1i * pc * (Y-q0(2)) / veps);
end


function [V, DV, D2V] = potential(Q1, Q2, quadraticScale, quarticScale)
    y1 = Q1 - 4;
    y2 = Q2 - 4;
    radiusSquared = y1.^2 + y2.^2;
    V = quadraticScale * radiusSquared + quarticScale * radiusSquared.^2;
    forceScale = 2 * quadraticScale + 4 * quarticScale * radiusSquared;
    DV = [forceScale .* y1, forceScale .* y2];
    H11 = forceScale + 8 * quarticScale * y1.^2;
    H22 = forceScale + 8 * quarticScale * y2.^2;
    H12 = 8 * quarticScale * y1 .* y2;
    D2V = [H11, H12, H12, H22];
end
