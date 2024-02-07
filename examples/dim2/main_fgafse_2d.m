%% Main function: compare FGA and TSSA for 2-dim fractional Schrodinger equation 
% with high-fequency wave

% add functions into file path
cd ../../
FGAFSE_startup();
cd ./examples/dim2

% ------------------------------------------------------------
% Various Scenarios for demonstration of error analysis
% ------------------------------------------------------------

% Option 1: for demonstration of alpha = 2, slope lies close to 1.0
% right_x = 2;
% final_time = 0.5;
% alpha = 2;
% vepsExp = [-6, -7, -8, -9];
% figureName = './L2err_alpha2.png';
% error_decay_2d(alpha, vepsExp, right_x, final_time, @initWave, @potential, figureName);


% Option 2: for demonstration of 5/3 < alpha < 2
% right_x = 2;
% final_time = 0.5;
% alpha = [1.7, 1.75, 1.80, 1.85, 1.90];
% vepsExp = [-6, -7, -8, -9];
% figureName = './L2err_alpha_5frac3_to_2.png';
% error_decay_2d(alpha, vepsExp, right_x, final_time, @initWave, @potential, figureName);


% Option 3: for demonstration of 1 < alpha < 2
right_x = 2;
final_time = 0.5;
alpha = [1.1, 1.3, 1.5, 1.7, 1.9];
vepsExp = [-6, -7, -8, -9];
figureName = './L2err_alpha_1to2.png';
error_decay_2d(alpha, vepsExp, right_x, final_time, @initWave, @potential, figureName);


% Option 4: for effect of exponential order k in delta = veps ^ k
% Note: default k=1, modify the value in odes_delta_2d() before running the test below.
% right_x = 2;
% final_time = 0.5;
% alpha = [1.1, 1.3, 1.5, 1.7, 1.9];
% vepsExp = [-6, -7, -8, -9];

% figureName = './L2err_alpha_1to2_k_2.png';
% figureName = './L2err_alpha_1to2_k_3/5.png';
% figureName = './L2err_alpha_1to2_delta_1e-8.png';

% error_decay_2d(alpha, vepsExp, right_x, final_time, @initWave, @potential, figureName);


% ------------------------------------------------------------
% Various choices of initial wavefunction and potential
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

    % V = 10;
    % DV_1 = zeros(size(Q1));
    % DV_2 = zeros(size(Q2));
    % DV = [DV_1, DV_2];
    % D2V = repmat([0, 0, 0, 0], size(Q1));
end