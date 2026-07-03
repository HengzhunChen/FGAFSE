%% Main function: compare FGA and TSSA for 2-dim fractional Schrodinger equation 
% with high-fequency wave

% add functions into file path
cd ../../
FGAFSE_startup();
cd ./examples/dim2

% ------------------------------------------------------------
% Demonstration of error analysis
% ------------------------------------------------------------

% Power k in delta = veps ^ k. The default used by FGA2d is k = 1.
deltaPower = 1;
figureName = './L2err_alpha_1to2.eps';
alpha = [1.1, 1.3, 1.5, 1.7, 1.9];

% deltaPower = 7/12;
% figureName = './L2err_alpha_1to2_k_7div12.eps';
% alpha = [1.75, 1.8, 1.85, 1.9, 1.95];

% Set false to recompute and overwrite TSSA data after changing its setup.
useCachedTSSA = true;

% % Demonstration of 1 < alpha < 2
right_x = 8;
final_time = 5;
vepsExp = [-5, -6, -7, -8];
error_decay_2d(alpha, vepsExp, right_x, final_time, @initWave, @potential, ...
               figureName, useCachedTSSA, deltaPower);


% ------------------------------------------------------------
% Initial wavefunction and potential
% ------------------------------------------------------------

function u0 = initWave(X, Y, veps)
% function to compute values of initial function
% inputs:
%       X, Y -- mesh samples
%       veps -- scaled Planck constant 
    
    r = sqrt( (X - 4).^2 + (Y - 4).^2 );
    u0 = sqrt(128 / pi) * exp(-64 * r.^2) .* exp(1i * (Y - 1) / veps);

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
    
    a = 0.1;
    b1 = 4; b2 = 4;
    V = a * ((Q1 - b1).^2 + (Q2 - b2).^2);
    DV = 2 * a * [Q1 - b1, Q2 - b2];
    D2V = repmat(2 * a * [1, 0, 0, 1], size(Q1));

    % V = 10;
    % DV_1 = zeros(size(Q1));
    % DV_2 = zeros(size(Q2));
    % DV = [DV_1, DV_2];
    % D2V = repmat([0, 0, 0, 0], size(Q1));
end
