function error_decay_1d(alpha, vepsExp, right_x, finalTime, initWave, potential, figName, useCachedTSSA, deltaPower)
% ERROR_DECAY_1D function to analyze the error decay rate of FGA for 1-dim
% fractional Schrodinger equation with high-frequency wave, use TSSA method as
% ground state.
%    Inputs:  
%        alpha     -- order of fractional operator
%        vepsExp   -- exponent of veps (scaled Planck constant), 
%                     veps = 2 ^ vepsExp
%        right_x   -- right endpoint of domain of x
%                     left endpoint of domain of x is 0
%        finalTime -- final time of evolution
%        initWave  -- function handle for initial wavefunction
%                     u0 = initWavefun(x, veps)
%        potential -- function handle of potential 
%                     [V, DV, D2V] = potential(Q)
%        figName   -- figure name for plot of L2 error decay curve
%        useCachedTSSA -- load saved TSSA results when available
%                         false recomputes and overwrites saved results
%        deltaPower -- power k in delta = veps ^ k
%                      default: deltaPower = 1
%
%    See also FGA1d, TSSA1d.

%  Copyright (c) 2024 Hengzhun Chen, Fudan University,
%                     Lihui Chai, Sun Yat-sen University.
%  This file is distributed under the terms of the MIT License.


% *****************************************************************************
%                        Error Analysis (+_+)
%
% (1) Convergence of TSSA
% 
% Since TSSA is used as ground truth, it should take small steps and at least 
% it should be convergent. At this time, it will not change a lot if we change 
% the time step to be finer, i.e.,
% 
%     || u(dt) - u(dt/2) || <= 1e-8.
%
% (2) Boundary condition
% 
% Note that fft() is used in TSSA method, we have assumed that it's periodic
% boundary condition, but we don't do this in FGA, thus it may not be periodic
% boundary. To avoid this, we should choose a suitable final time to restrict
% the wave not to go outside the space domain.
%
% (3) Error in FGA procedure and choices of parameters in FGA
%
%     err_comp = err_FGA + err_odes + err_init
%
% If we want to demonstrate relation between veps and err_FGA, we have to make
% err_comp = err_FGA, i.e., let err_odes and err_init be small enough to be
% ignored. For instance,
% 
%     err_odes, err_init <= 1e-16, 
%
% which is the machine error; for another, at least, 
%
%     err_odes = O(err_FGA) and err_init = O(err_FGA)
%
% We should choose proper parameters(mesh strategy), especially their relation
% with respect to veps.
% For err_odes, RK4 has fourth-order global accuracy, while the symplectic
% Stormer-Verlet solver has second-order global accuracy.
%
%     err_odes = O(dt^4) for RK4, O(dt^2) for Stormer-Verlet
%
% thus, we can take dt = 1e-3 or even dt = 1e-2.
% For err_init, since we set dx = veps, it will get smaller whenever veps get
% smaller, at least, we ask other mesh strategy such that
%
%     err_init = O(err_comp)
%
% *****************************************************************************

if nargin < 8
    useCachedTSSA = false;
end
if nargin < 9 || isempty(deltaPower)
    deltaPower = 1;
end

% FGA ODE solver: 'symplectic' or 'rk4'.
% solver = 'symplectic';
solver = 'rk4';

nalpha = length(alpha);  % number of alpha to test
nveps = length(vepsExp);  % number of veps to test

err_L2 = zeros(nalpha, nveps);
time_FGA = zeros(nalpha, nveps);
time_TSSA = zeros(nalpha, nveps);

folder = './figures/error_decay';
if ~exist(folder, 'file')
    mkdir(folder);
end
tssaFolder = './data/tssa_1d';
if ~exist(tssaFolder, 'dir')
    mkdir(tssaFolder);
end
old_png = folder + "/alpha*veps*.png";
old_eps = folder + "/alpha*veps*.eps";
delete(old_png);  % Note: clean old figures
delete(old_eps);  % Note: clean old figures

for i = 1 : nalpha
    for j = 1 : nveps
        veps = 2 ^ vepsExp(j);
        dx = 2 ^ vepsExp(j);
        fprintf("alpha = %f, veps = %e\n", alpha(i), veps);

        % -------------------------------------------------
        % Solver computation
        % -------------------------------------------------

        % FGA method, display error of initial decompsition
        timer_FGA = tic;
        delta = veps ^ deltaPower;
        [w, x_w] = FGA1d(alpha(i), vepsExp(j), finalTime, right_x, ...
            initWave, potential, delta, solver);
        time_FGA(i, j) = toc(timer_FGA);
        fprintf('FGA1d time: %f seconds\n', time_FGA(i, j));

        % TSSA method, as ground truth
        % dt = veps;
        dt = veps ^ 2;
        cacheFile = fullfile(tssaFolder, sprintf( ...
            'alpha_%g_vepsExp_%d_T_%g_R_%g.mat', ...
            alpha(i), vepsExp(j), finalTime, right_x));
        timer_TSSA = tic;
        if useCachedTSSA && exist(cacheFile, 'file')
            cached = load(cacheFile, "u", "x_u");
            u = cached.u;
            x_u = cached.x_u;
            fprintf('Loaded TSSA1d result: %s\n', cacheFile);
        else
            [u, x_u] = TSSA1d(alpha(i), vepsExp(j), finalTime, right_x, dt, dx, initWave, potential);
            save(cacheFile, "u", "x_u");
            fprintf('Saved TSSA1d result: %s\n', cacheFile);
        end
        time_TSSA(i, j) = toc(timer_TSSA);
        fprintf('TSSA1d result time: %f seconds\n', time_TSSA(i, j));

        % test for convergence of TSSA method
        % dt = dt / 2;
        % [u1, ~] = TSSA1d(alpha(i), vepsExp(j), finalTime, right_x, dt, dx, initWave, potential);
        % dist_TSSA = norm(u-u1) / norm(u);
        % fprintf('distance between two time steps of TSSA is %e\n', dist_TSSA);

        % -----------------------------------------------
        % Visualization
        % -----------------------------------------------

        % ----  Code for plot figures used in paper ----
        % figure;
        % xidx = x_w < 1.8 & x_w > 0.8; 
        % hold on
        % box on
        % plot(x_u(xidx), real(u(xidx)), 'b-', 'LineWidth', 2.5);
        % plot(x_w(xidx), real(w(xidx)), 'r--', 'LineWidth', 2.5);
        % xlabel('x', 'FontSize', 16, 'FontWeight', 'bold')
        % ylabel('Re \psi', 'FontSize', 16, 'FontWeight', 'bold')
        % ylim([-6, 6])
        % set(gca, 'Fontsize', 15, 'FontWeight', 'bold')
        % hold off
        % filename = folder + "/alpha_" + num2str(i) + "_veps_" + num2str(j);
        % saveas(gcf, filename, 'png');
        % saveas(gcf, filename, 'epsc');

        % ---- Code for further solution comparison ----
        figure;
        subplot(2, 2, 1)
        hold on
        plot(x_w, real(w), '-');
        plot(x_u, real(u), '-.');
        hold off
        title('real part')

        subplot(2, 2, 2)
        hold on
        plot(x_w, imag(w), '-');
        plot(x_u, imag(u), '-.');
        hold off
        title('imag part')

        subplot(2, 2, 3)
        hold on
        plot(x_w, abs(w).^2, '-');
        plot(x_u, abs(u).^2, '-.');
        hold off
        title('position density')

        subplot(2, 2, 4)
        hold on    
        plot(x_w, veps * imag( conj(w) .* gradient(w) ), '-');
        plot(x_u, veps * imag( conj(u) .* gradient(u) ), '-.');
        hold off
        title('current density')

        legend('FGA','TSSA', 'Orientation', 'horizontal', 'Location', [0.52 0.03  0  0])
        sgtitle(['alpha = ', num2str(alpha(i)), ', t = ', num2str(finalTime), ...
                ', varepsilon = ', num2str(veps)]);
        filename = folder + "/alpha_" + num2str(i) + "_veps_" + num2str(j);
        saveas(gcf, filename, 'png');

        % ------------------------------------------------
        % Error calculation
        % ------------------------------------------------
        err_L2(i, j) = sqrt( sum( abs(u - w).^2 ) * dx );
        fprintf('L2 distance between FGA and TSSA: %e\n', err_L2(i, j));
    end
end

% ------------------------------------------
% Print result of errors
% ------------------------------------------
fprintf("\nTable of L2 error\n");
for i = 1: nalpha
    for j = 1 : nveps
        fprintf('%.2e ', err_L2(i, j));
    end
    fprintf('\n')
end

fprintf("\nTable of FGA1d time in seconds\n");
for i = 1: nalpha
    for j = 1 : nveps
        fprintf('%.2e ', time_FGA(i, j));
    end
    fprintf('\n')
end

fprintf("\nTable of TSSA1d result time in seconds\n");
for i = 1: nalpha
    for j = 1 : nveps
        fprintf('%.2e ', time_TSSA(i, j));
    end
    fprintf('\n')
end

save("L2_error.mat", "err_L2", "time_FGA", "time_TSSA");

% --------------------------------------
% Plot error decay curves
% --------------------------------------

% Fix alpha to plot the relation between log(error) and log(veps),
% hopefully, error = O(veps ^ r), find out r
mylinestyle = ["-*", "-o", "-^", "-square", "-diamond"];
figure;
set(gca, 'Fontsize', 14)
hold on
box on 
grid on

p1 = zeros(nalpha, 2);  % coefficents of ployfit
leg_str = cell(nalpha, 1);
for i = 1 : nalpha
    plot(-vepsExp, log2(err_L2(i, :)), mylinestyle(i), 'Linewidth', 1);
    p1(i, :) = polyfit(vepsExp, log2(err_L2(i, :)), 1);
    leg_str{i} = ['$\alpha=$', num2str(alpha(i)), ', slope: ', num2str(p1(i, 1))];
end

txt_x = xlabel("$-\log_2(\varepsilon)$");
txt_y = ylabel("$\log_2$($L^2$ error)");
set(txt_x, "Interpreter", "latex")
set(txt_y, "Interpreter", "latex")
legend(leg_str, "Interpreter", "latex");
% title(['L^2 error with final time T = ', num2str(finalTime)]);
hold off

if isempty(figName)
    figName = './L2_err_veps.png';
end
if contains(figName, "eps")
    saveas(gcf, figName, 'epsc');
    figName = replace(figName, "eps", "png");
end
saveas(gcf, figName, 'png');

fprintf('    slope    intercept  (L2 error)\n');
disp(p1);

end
