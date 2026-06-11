%% Test for functions of 1-dim FGA components work correctly

% add functions into file path
cd ../
FGAFSE_startup();
cd ./test

right_x = 3;
final_time = 2;

vepsExp = -9;
veps = 2 ^ vepsExp;
alpha = 1.2;
dx = veps;

initWave = @initWavefun;
potential = @potentialfun;

[w, x_w] = FGA1d(alpha, vepsExp, final_time, right_x, initWave, potential);

dt = veps^2;
[u, x_u] = TSSA1d(alpha, vepsExp, final_time, right_x, dt, dx, initWave, potential);

% Visualization
folder = './figures';
if ~exist(folder, 'file')
    mkdir(folder);
end
figure;

subplot(2, 2, 1);
hold on
plot(x_w, real(w), '-');
plot(x_u, real(u), '-.');
hold off
title('Subplot 1: real(u)')

s2 = subplot(2, 2, 2);
hold on
plot(x_w, imag(w), '-');
plot(x_u, imag(u), '-.');
hold off
title('Subplot 2: imag(u)')

s3 = subplot(2,2,3);
hold on
plot(x_w, abs(w) .^ 2, '-');
plot(x_u, abs(u) .^ 2, '-.');
hold off
title('Subplot 3: |u|^2')

s4 = subplot(2,2,4);
hold on
% plot(x_w, imag( gradient(w .* conj(w)) ), '-');
% plot(x_u, imag( gradient(u .* conj(u)) ), '-.');
plot(x_w, veps * imag( conj(w) .* gradient(w) ), '-');
plot(x_u, veps * imag( conj(u) .* gradient(u) ), '-.');
hold off
title('Subplot 4: current density')

legend('FGA','TSSA', 'Orientation', 'horizontal', 'Location', [0.52 0.03  0  0])

errL2 = sqrt( sum( abs(w - u).^2 ) * dx );
fprintf('L2 error between FGA and TSSA: %e\n', errL2);

sgtitle(['\alpha = ', num2str(alpha), ', t = ', num2str(final_time), ... 
         ', L2 error = ', num2str(errL2)]);  

saveas(gcf, './figures/test1d.png', 'png');

% compute the L2 norm of initial wavefunction
x_init = linspace(0, right_x, 1000);
u0_init = initWave(x_init, veps);
L2_norm_u0 = sqrt(sum(abs(u0_init).^2) * (x_init(2) - x_init(1)));
fprintf('L2 norm of initial wavefunction: %e\n', L2_norm_u0);

% compute the L2 norm of the solution
L2_norm_w = sqrt(sum(abs(w).^2) * dx);
L2_norm_u = sqrt(sum(abs(u).^2) * dx);
fprintf('L2 norm of FGA solution: %e\n', L2_norm_w);
fprintf('L2 norm of TSSA solution: %e\n', L2_norm_u);


% Test for counting zero points in trajectory P(t)
count_zero_P_1d(alpha, final_time, vepsExp, right_x, @initWavefun, @potentialfun);


% ------------------------------------------------------------

function u0 = initWavefun(x, veps)
    beta = 1;
    u0 = exp(-64 * (x - 0.5).^2) / sqrt(pi / 64) .* exp(1i / veps * beta * x * 1);
end


function [V, DV, D2V] = potentialfun(Q)
    a = 1;
    b = 1.5;
    V = a * (Q - b).^2;
    DV = 2 * a * (Q - b);
    D2V = 2 * a;
end
