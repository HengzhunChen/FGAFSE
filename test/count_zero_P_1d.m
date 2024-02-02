function cnt1 = count_zero_P_1d(alpha, finalTime, vepsExp, right_x, initWave, potential)
% COUNT_ZERO_P_1D Count the number of zero points in the trajectory of P(t) in 1 dimension.
%    Inputs:
%        alpha     -- order of fractional operator
%        finalTime -- final time of evolution
%        vepsExp   -- exponent of veps (scaled Planck constant), 
%                     veps = 2 ^ vepsExp
%        right_x   -- right endpoint of domain of x
%                     left endpoint of domian of x is 0
%        initWave  -- function handle for initial wavefunction
%                     u0 = initWavefun(x, veps)
%        potential -- function handle of potential 
%                     [V, DV, D2V] = potential(Q)
%    Outputs:
%       cnt1 -- number of zero points in the trajectory of P(t) 
%
%    See also initial_deomposition_1d, odes_delta_1d, wave_reconstruction_1d.


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

dt = 1e-2;
nt = floor( finalTime / dt );

% Initialization
x = 0 : dx : right_x;  % mesh on axis x, left endpoint is 0
x = x(1 : end-1)';  % shape: (nx, 1)
u0 = initWave(x, veps);

odes = @odes_delta_1d;

% Time evolution
[A0, S0, Q0, P0, nGB] = initial_decomposition_1d(u0, veps, dy, ny, kernelSize, nydq);
DzQ0 = ones(size(Q0));
DzP0 = -1i * ones(size(P0));
Q = zeros(nGB, nt);
P = zeros(nGB, nt);
for tt = 1 : nt
    [A0, S0, Q0, P0, DzQ0, DzP0] = time_evolution(A0, S0, Q0, P0, DzQ0, DzP0, dt, 1, alpha, odes, potential, veps);
    Q(:, tt) = Q0;
    P(:, tt) = P0;
end

cnt1 = 0;  % count for pairs (q,p) such that P(t,q,p) go through zero points
cnt2 = 0;  % count for singular points of the Hamiltonian

% Visualization of FGA quantities

% folder = './figures/QP_demo_1d';
% if ~exist(folder, 'file')
%     mkdir(folder);
% end
% t = 0 : dt : finalTime;
% t = t(1 : end-1);

for i = 1 : nGB
    P_qp = P(i, :);
    Q_qp = Q(i, :);

    [~, DV, ~] = potential(Q_qp);

    if (max(P_qp) > 1e-12) && (min(P_qp) < -1e-12)
        cnt1 = cnt1 + 1;

        % plot those Q(t,q,p), P(t,q,p) curves such that P(t) go through zero points
        % figure;
        % subplot(1, 2, 1);
        % plot(t, P_qp, '-');
        % title('P')
        % subplot(1, 2, 2);
        % plot(t, Q_qp, '-');
        % title('Q')
        % saveas(gcf, "./figures/QP_demo_1d/qp_" + num2str(i) + "_.png", 'png');

        PQ = abs(P_qp) + abs(DV);
        if min(PQ) < 1e-4
            cnt2 = cnt2 + 1;
        end
    end
end
fprintf('During time interval [0, %.2f],\n', finalTime)
fprintf('there are %d out of %d initial value pairs (q,p) such that P(t,q,p) go through zero points,\n', cnt1, nGB);
fprintf('there are %d singular points\n', cnt2);

end