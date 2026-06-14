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
figureName = './L2err_alpha_1to2.png';
% deltaPower = 7/12;
% figureName = './L2err_alpha_1to2_k_7div12.png';

% Set false to recompute and overwrite TSSA data after changing its setup.
useCachedTSSA = false;

% % Demonstration of 1 < alpha < 2
right_x = 2;
final_time = 2;
alpha = [1.1, 1.3, 1.5, 1.7, 1.9];
vepsExp = [-6, -7, -8, -9];
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
    
    V = ((Q1 - 0.5).^2 + (Q2 - 0.5).^2) / 2;
    DV_1 = Q1 - 0.5;
    DV_2 = Q2 - 0.5;    
    DV = [DV_1, DV_2];
    D2V = repmat([1, 0, 0, 1], size(Q1));

    % V = 10;
    % DV_1 = zeros(size(Q1));
    % DV_2 = zeros(size(Q2));
    % DV = [DV_1, DV_2];
    % D2V = repmat([0, 0, 0, 0], size(Q1));
end
