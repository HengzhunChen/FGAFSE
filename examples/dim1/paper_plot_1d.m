%% Plot figures used in the paper

% add functions into file path
cd ../../
FGAFSE_startup();
cd ./examples/dim1

% NOTE: remember to comment the visualization code in error_decay_1d() to match the requirement.

% Solution Comparison
% right_x = 2;
% final_time = 0.25;
% alpha = [1.3, 1.7];
% vepsExp = [-6, -7, -8];
% figureName = './L2err_alpha_1to2.png';
% error_decay_1d(alpha, vepsExp, right_x, final_time, @initWave, @potential, figureName);

% Error decay
right_x = 2;
final_time = 0.25;
alpha = [1.1, 1.3, 1.5, 1.7, 1.9];
vepsExp = [-6, -7, -8, -9, -10];
figureName = './L2err_alpha_1to2.png';
% figureName = './L2err_alpha_1to2_k_1.eps';
% figureName = './L2err_alpha_1to2_k_6frac11.eps';
error_decay_1d(alpha, vepsExp, right_x, final_time, @initWave, @potential, figureName);
for alphaIdx = alpha
    for expIndex = vepsExp
        % check whether the trajectory of P go through its zero points
        fprintf("alpha = %f, vepsExp = %d\n", alphaIdx, expIndex);
        count_zero_P_1d(alphaIdx, final_time, expIndex, right_x, @initWave, @potential);
    end
end


% ------------------------------------------------------------
% Initial wavefunction and potential
% ------------------------------------------------------------

function u0 = initWave(x, veps)
% function to compute values of initial function
% inputs:
%       veps -- scaled Planck constant 
    beta = 1;
    u0 = exp(-64 * (x - 1).^2) / sqrt(pi / 64) .* exp(1i / veps * beta * x * 1);
end


function [V, DV, D2V] = potential(Q)
% function to compute values and derivatives of potential function 
% input:
%        Q -- position where the potential to be evaluated
% outputs:
%        V = V(Q)    potential value 
%        DV = V'(Q)   1st derivative of V
%        D2V = V''(Q)  2nd derivative of V
    V = 1 + cos(pi * Q);
    DV = -pi * sin(pi * Q);
    D2V = -(pi)^2 * cos(pi * Q);
end
