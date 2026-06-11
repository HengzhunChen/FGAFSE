%% Main function: compare FGA and TSSA for 1-dim fractional Schrodinger equation 
%% with high-frequency wave in semi-classical regime

% add functions into file path
cd ../../
FGAFSE_startup();
cd ./examples/dim1

% ------------------------------------------------------------
% Demonstration of error analysis
% ------------------------------------------------------------

% % Power k in delta = veps ^ k. The default used by FGA1d is k = 1.
deltaPower = 1;
figureName = './L2err_alpha_1to2.png';
% deltaPower = 7/12;
% figureName = './L2err_alpha_1to2_k_7div12.png';

% Set false to recompute and overwrite TSSA data after changing its setup.
useCachedTSSA = true;

% % Demonstration of 1 < alpha < 2
right_x = 3;
final_time = 2;
alpha = [1.1, 1.3, 1.5, 1.7, 1.9];
vepsExp = [-6, -7, -8, -9, -10];
error_decay_1d(alpha, vepsExp, right_x, final_time, @initWave, @potential, ...
               figureName, useCachedTSSA, deltaPower);

% % Check whether the trajectories of P pass through zero.
% for alphaIdx = alpha
%     for expIndex = vepsExp
%         fprintf("alpha = %f, vepsExp = %d\n", alphaIdx, expIndex);
%         count_zero_P_1d(alphaIdx, final_time, expIndex, right_x, @initWave, @potential);
%     end
% end


% ------------------------------------------------------------
% Initial wavefunction and potential
% ------------------------------------------------------------

function u0 = initWave(x, veps)
% function to compute values of initial function
% inputs:
%       veps -- scaled Planck constant 
    beta = 1;
    u0 = (128 / pi)^(1/4) * exp(-64 * (x - 1.5).^2) ...
        .* exp(1i / veps * beta * x);
end


function [V, DV, D2V] = potential(Q)
% function to compute values and derivatives of potential function 
% input:
%        Q -- position where the potential to be evaluated
% outputs:
%        V = V(Q)    potential value 
%        DV = V'(Q)   1st derivative of V
%        D2V = V''(Q)  2nd derivative of V
    a = 1;
    b = 1.5;
    V = a * (Q - b).^2;
    DV = 2 * a * (Q - b);
    D2V = 2 * a;

    % This potential will show relation between initial error and veps 
    % since both FGA and TSSA are exact.
    % V = 10;
    % DV = 0;
    % D2V = 0;
end
