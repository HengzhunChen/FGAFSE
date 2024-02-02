function error_decay_1d(alpha, vepsExp, right_x, finalTime, initWave, potential, figName)
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
% For err_odes, since we use 4-th order Runge-Kutta method, 
%
%     err_odes = O(dt^5)
%
% thus, we can take dt = 1e-3 or even dt = 1e-2.
% For err_init, since we set dx = veps, it will get smaller whenever veps get
% smaller, at least, we ask other mesh strategy such that
%
%     err_init = O(err_comp)
%
% *****************************************************************************

nalpha = length(alpha);  % number of alpha to test
nveps = length(vepsExp);  % number of veps to test

err_L2 = zeros(nalpha, nveps);

folder = './figures/error_decay';
if ~exist(folder, 'file')
    mkdir(folder);
end

for i = 1 : nalpha
    for j = 1 : nveps
        veps = 2 ^ vepsExp(j);
        dx = 2 ^ vepsExp(j);
        fprintf("alpha = %f, veps = %e\n", alpha(i), veps);

        % -------------------------------------------------
        % Solver computation
        % -------------------------------------------------

        % FGA method, display error of initial decompsition
        [w, x_w] = FGA1d(alpha(i), vepsExp(j), finalTime, right_x, initWave, potential);

        % TSSA method, as ground truth
        % dt = veps;
        dt = veps ^ 2;
        [u, x_u] = TSSA1d(alpha(i), vepsExp(j), finalTime, right_x, dt, dx, initWave, potential);

        % test for convergence of TSSA method
        dt = dt / 2;
        [u1, ~] = TSSA1d(alpha(i), vepsExp(j), finalTime, right_x, dt, dx, initWave, potential);
        dist_TSSA = norm(u-u1) / norm(u);
        fprintf('distance between two time steps of TSSA is %e\n', dist_TSSA);

        % -----------------------------------------------
        % Visualization
        % -----------------------------------------------
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

        saveas(gcf, ['./figures/error_decay/', 'alpha_', num2str(i), 'veps_', num2str(j), '.png']);

        % ------------------------------------------------
        % Error calculation
        % ------------------------------------------------
        err_L2(i, j) = sqrt( sum( abs(u - w).^2 ) * dx );
        fprintf('L2 distance between FGA and TSSA: %e\n', err_L2(i, j));
    end
end

% ------------------------------------------
% Write result to a table
% ------------------------------------------
colname = strings(1, nveps);
for j = 1 : nveps
    colname(1, j) = append("veps=1/", num2str(2^(-vepsExp(j))));
end
rowname = strings(nalpha, 1);
for i = 1 : nalpha
    rowname(i, 1) = append("alpha=", num2str(alpha(i)));
end
L2Table = ["L2 error", colname; rowname, err_L2];
fprintf("\nTable of L2 error\n");
disp(L2Table);
% writematrix(L2Table, 'error.xlsx');

% --------------------------------------
% Plot error decay curves
% --------------------------------------

% Fix alpha to plot the relation between log(error) and log(veps),
% hopefully, error = O(veps ^ r), find out r
figure;
hold on
p1 = zeros(nalpha, 2);  % coefficents of ployfit
leg_str = cell(nalpha, 1);
for i = 1 : nalpha    
    plot(-vepsExp, log2(err_L2(i, :)), '-*');
    p1(i, :) = polyfit(vepsExp, log2(err_L2(i, :)), 1);
    leg_str{i} = ['alpha=', num2str(alpha(i)), ', slope: ', num2str(p1(i, 1))];
    xlabel("-log_2(veps)");
    ylabel("log_2(L2 error)");
end
legend(leg_str);
title(['L2 error with final time T = ', num2str(finalTime)]);
if isempty(figName)
    figName = './L2_err_veps.png';
end
saveas(gcf, figName);
hold off

fprintf('    slope    intercept  (L2 error)\n');
disp(p1);

end