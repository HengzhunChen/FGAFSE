%% Plot figures used in the paper

% add functions into file path
cd ../../
FGAFSE_startup();
cd ./examples/dim2

% Solution Comparison
right_x = 2;
final_time = 0.25;
alpha = 1.5;
vepsExp = [-6, -7, -8];
figureName = './L2err_alpha_1to2.png';
error_decay_2d(alpha, vepsExp, right_x, final_time, @initWave, @potential, figureName);

% Error decay
% right_x = 2;
% final_time = 0.25;
% alpha = [1.1, 1.3, 1.5, 1.7, 1.9];
% vepsExp = [-6, -7, -8, -9];
% % figureName = './L2err_alpha_1to2.png';
% figureName = './L2err_alpha_1to2.eps';
% error_decay_2d(alpha, vepsExp, right_x, final_time, @initWave, @potential, figureName);


% ------------------------------------------------------------
% Initial wavefunction and potential
% ------------------------------------------------------------

function u0 = initWave(X, Y, veps)
% function to compute values of initial function
% inputs:
%       X, Y -- mesh samples
%       veps -- scaled Planck constant     
    r = sqrt( (X - 1).^2 + (Y - 1).^2 );
    u0 = exp( -(r.^2) * 64 ) / (pi / 64) .* exp( 1i * (Y - 1) / veps);
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
    V = ((Q1 - 1) .^ 2 + (Q2 - 1) .^ 2) / 2;
    DV_1 = Q1 - 1;
    DV_2 = Q2 - 1;    
    DV = [DV_1, DV_2];
    D2V = repmat([1, 0, 0, 1], size(Q1));
end