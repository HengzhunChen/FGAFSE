function check_zero_P_2d(alpha, finalTime, vepsExp, right_x, delta, initWave, potential)
% CHEcK_ZERO_P_2D Check the existence of zero points in the trajectory of |P(t)| in 2 dimension.
%    Inputs:
%        alpha     -- order of fractional operator
%        finalTime -- final time of evolution
%        vepsExp   -- exponent of veps (scaled Planck constant), 
%                     veps = 2 ^ vepsExp
%        right_x   -- right endpoint of domain of x
%                     left endpoint of domian of x is 0
%        delta     -- parameter for modified Laplacian operator
%        initWave  -- function handle for initial wavefunction
%                     u0 = initWave(X, Y, veps)
%        potential -- function handle of potential 
%                     [V, DV, D2V] = potential(Q1, Q2)
%
%    See also initial_deomposition_2d.


veps = 2 ^ vepsExp;  % scaled Planck constant

% Setup mesh grid
dx = veps;
nx = floor( (right_x - 0) / dx);
dy = dx;  % use the same mesh size for x and y
ny = nx;
% number of y grid included in each stepsize of q, nydq := ny / nq, dq := dy * nqdq
nydq = floor( 2^(-vepsExp / 2) / 2 );
% number of points included in a Gaussian kernel
kernelSize = floor( 2^(-vepsExp / 2) ) * 2^3;

% Initialization
x = linspace(0, right_x, nx+1);  % shape: (1, nx+1)
x = x(1 : end-1)';  % shape: (nx, 1)
xx = x * ones(1, nx);
u0 = initWave(xx, xx', veps);
[~, ~, Qmesh, Pmesh, ~] = initial_decomposition_2d(u0, veps, dy, ny, kernelSize, nydq);

P1 = Pmesh(:, 1); P2 = Pmesh(:, 2);
Q1 = Qmesh(:, 1); Q2 = Qmesh(:, 2);

P1min = min(P1); P1max = max(P1);
P2min = min(P2); P2max = max(P2);
Q1min = min(Q1); Q1max = max(Q1);
Q2min = min(Q2); Q2max = max(Q2);

[~, DV, ~] = potential(Q1, Q1);
DV_norm = sqrt( DV(:, 1).^2 + DV(:, 2).^2 );
index = DV_norm > 1e-8;
Q1 = Q1(index);
Q2 = Q2(index);

% Setup end condition for inverse ODEs
Qi = [Q1, Q2];
Pi = zeros(size(Qi));

dt = 1e-2;
halfTime = finalTime / 2;
Nt = floor(halfTime / dt + 1e-6);

% Main loop
for tt = 1 : Nt
    % find slope vector k1
    [kQ1, kP1] = inverse_odes2d(Qi, Pi, alpha, potential, delta);
    
    % find slope vector k2
    Q = Qi + kQ1 * dt / 2;
    P = Pi + kP1 * dt / 2;
    [kQ2, kP2] = inverse_odes2d(Q, P, alpha, potential, delta);
    
    % find slope vector k3       
    Q = Qi + kQ2 * dt / 2;
    P = Pi + kP2 * dt / 2;
    [kQ3, kP3] = inverse_odes2d(Q, P, alpha, potential, delta);
    
    % find slope vector k4
    Q = Qi + kQ3 * dt;
    P = Pi + kP3 * dt;
    [kQ4, kP4] = inverse_odes2d(Q, P, alpha, potential, delta);

    % evolution
    Qi = Qi + (kQ1 + 2 * kQ2 + 2 * kQ3 + kQ4) * dt / 6;
    Pi = Pi + (kP1 + 2 * kP2 + 2 * kP3 + kP4) * dt / 6;
end

Pi1 = Pi(:, 1); Pi2 = Pi(:, 2);
Qi1 = Qi(:, 1); Qi2 = Qi(:, 2);

idx = P1min < Pi1 & Pi1 < P1max & P2min < Pi1 & Pi2 < P2max & ...
      Q1min < Qi1 & Qi1 < Q1max & Q2min < Qi2 & Qi2 < Q2max;

% Count for those trajectories go through |P|=0 and Q inside grid domain at half
% time while their initial conditions are also inside the grid domain.
cnt = length(idx);

if cnt > 10
    fprintf('There are %d pair target initial conditions satisfying domain constraint.\n', cnt);
    fprintf('There exists some trajectories such that P(t,q,p) go through zero points during evolution.\n');
end

end


% ---------------------------------------------------------------------
% Subfunctions for the righthand side of ODEs in inverse time direction
% ---------------------------------------------------------------------

function [DQ, DP] = inverse_odes2d(Q, P, alpha, potential, delta)    
    [~, DV, ~] = potential(Q(:, 1), Q(:, 2));
    
    if alpha == 2
        DQ = -P;
        DP = DV;        
    else    
        Pdelta = sum(P.^2, 2) + delta^2;
        DQ = -Pdelta .^ ((alpha-2)/2) .* P;
        DP = DV;
    end
end
