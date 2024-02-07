%% Main function: compare FGA and TSSA for 1-dim fractional Schrodinger equation 
% with high-frequency wave

% add functions into file path
cd ../../
FGAFSE_startup();
cd ./examples/dim1

% ------------------------------------------------------------
% Various Scenarios for demonstration of error analysis
% ------------------------------------------------------------

% Option 1: for demonstration of alpha = 2, slope lies close to 1.0
% right_x = 2;
% final_time = 0.5;
% alpha = 2;
% vepsExp = [-6, -7, -8, -9, -10];
% figureName = './L2err_alpha2.png';
% error_decay_1d(alpha, vepsExp, right_x,  final_time, @initWave, @potential, figureName);
% for expIndex = vepsExp
%     % check whether the trajectory of P go through its zero points
%     fprintf("alpha = %f, vepsExp = %d\n", alpha, expIndex);
%     count_zero_P_1d(alpha, final_time, expIndex, right_x, @initWave, @potential);
% end

% Option 2: for demonstration of 11/6 < alpha < 2
% right_x = 2;
% final_time = 0.5;
% alpha = [1.875, 1.9, 1.925, 1.95, 1.975];
% vepsExp = [-6, -7, -8, -9, -10];
% figureName = './L2err_alpha_11frac6_to_2.png';
% error_decay_1d(alpha, vepsExp, right_x, final_time, @initWave, @potential, figureName);
% for alphaIdx = alpha
%     for expIndex = vepsExp
%         % check whether the trajectory of P go through its zero points
%         fprintf("alpha = %f, vepsExp = %d\n", alphaIdx, expIndex);
%         count_zero_P_1d(alphaIdx, final_time, expIndex, right_x, @initWave, @potential);
%     end
% end

% Option 3: for demonstration of 1 < alpha < 2
right_x = 2;
final_time = 0.5;
alpha = [1.1, 1.3, 1.5, 1.7, 1.9];
vepsExp = [-6, -7, -8, -9, -10];
figureName = './L2err_alpha_1to2.png';
error_decay_1d(alpha, vepsExp, right_x, final_time, @initWave, @potential, figureName);
for alphaIdx = alpha
    for expIndex = vepsExp
        % check whether the trajectory of P go through its zero points
        fprintf("alpha = %f, vepsExp = %d\n", alphaIdx, expIndex);
        count_zero_P_1d(alphaIdx, final_time, expIndex, right_x, @initWave, @potential);
    end
end

% Option 4: for effect of exponential order k in delta = veps ^ k
% Note: default k=1, modify the value in odes_delta_1d() before running the test below.
% right_x = 2;
% final_time = 0.5;
% alpha = [1.1, 1.3, 1.5, 1.7, 1.9];
% vepsExp = [-6, -7, -8, -9, -10];

% figureName = './L2err_alpha_1to2_k_2.png';
% figureName = './L2err_alpha_1to2_k_11/6.png';
% figureName = './L2err_alpha_1to2_delta_1e-8.png';

% error_decay_1d(alpha, vepsExp, right_x, final_time, @initWave, @potential, figureName);
% for alphaIdx = alpha
%     for expIndex = vepsExp
%         % check whether the trajectory of P go through its zero points
%         fprintf("alpha = %f, vepsExp = %d\n", alphaIdx, expIndex);
%         count_zero_P_1d(alphaIdx, final_time, expIndex, right_x, @initWave, @potential);
%     end
% end


% ------------------------------------------------------------
% Various choices of initial wavefunction and potential
% ------------------------------------------------------------

function u0 = initWave(x, veps)
% function to compute values of initial function
% inputs:
%       veps -- scaled Planck constant 

    beta = 1;
    u0 = exp(-64 * (x - 0.5).^2) / sqrt(pi / 64) .* exp(1i / veps * beta * x * 1);

    % n0 = (exp( -25 * (x - 0.5).^2 )) .^ 2;
    % S0 = x + 1;    
    % u0 = sqrt(n0) .* exp( 1i * S0 / veps );
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

    % V = exp(Q);
    % DV = exp(Q);
    % D2V = exp(Q);

    % V = 0.5 * ( (Q-1.0).^2 );
    % DV = Q-1.0;
    % D2V = 1.0;

    % This potential will show relation between initial error and veps 
    % since both FGA and TSSA are exact.
    % V = 10;
    % DV = 0;
    % D2V = 0;
end
