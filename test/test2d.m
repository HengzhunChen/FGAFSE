%% Test for functions of 2-dim FGA components work correctly

% parameters setting
right_x = 2;
final_time = 0.5;

vepsExp = -6;
veps = 2 ^ vepsExp;
alpha = 1.8;

% FGA method
[w, xx_w] = FGA2d(alpha, vepsExp, final_time, right_x, @initWave, @potential);

% time splitting method
dt = veps ^ 2;
[u, xx_u] = TSSA2d(alpha, vepsExp, final_time, right_x, dt, @initWave, @potential);


% Visualization
folder = './figures';
if ~exist(folder, 'file')
    mkdir(folder);
end
figure;
subplot(1, 2, 1)
hold on
pcolor(xx_u, xx_u', abs(u) .^ 2);
shading interp;
hold off
title('Subplot 1: solution of TSSA')

subplot(1, 2, 2)
hold on
pcolor(xx_w, xx_w', abs(w) .^ 2);
shading interp;
hold off
title('Subplot 2: solution of FGA')

% error computation
dx = veps;
err_L2 = sqrt( sum(sum( abs(u - w).^2 )) * dx * dx );
disp(['L2 error: ', num2str(err_L2)]);
err_L2 = sqrt( sum( sum( abs(u - w).^2 )) / sum(sum( abs(u).^2 )) );
disp(['normalized L2 error: ', num2str(err_L2)]);

saveas(gcf, './figures/test2d.png', 'png');

% Test for checking zero points in trajectory P(t)
delta = veps;
check_zero_P_2d(alpha, final_time, vepsExp, right_x, delta, @initWave, @potential);


% ------------------------------------------------------------

function u0 = initWave(X, Y, veps)
% function to compute values of initial function
% inputs:
%       X, Y -- mesh samples
%       veps -- scaled Planck constant 
    r = sqrt( (X - 0.5).^2 + (Y - 0.5).^2 );
    u0 = exp( -(r.^2) * 64 ) / (pi / 64) .* exp( 1i * (Y - 0.5) / veps);
end


function [V, DV, D2V] = potential(Q1, Q2)
% function to compute values and derivatives of potential function 
% input:
%        Q1, Q2 -- independent variables, can be vector or matrix 
% outputs:
%        V   -- V(Q1, Q2), potential value 
%        DV  -- [DV_1, DV_2]
%               1st partial derivatives of V with respect to q1, q2
%        D2V -- [DV_11, DV_12, DV_21, DV_22]
%               2nd partial derivatives of V w.r.t q1, q1
    V = ((Q1 - 0.5) .^ 2 + (Q2 - 0.5) .^ 2) / 2;
    DV_1 = Q1 - 0.5;
    DV_2 = Q2 - 0.5;    
    DV = [DV_1, DV_2];
    D2V = repmat([1, 0, 0, 1], size(Q1));
    
    % c = 10;  % scaling coefficient
    % V = V * c;
    % DV = DV * c;
    % D2V = D2V * c;
end