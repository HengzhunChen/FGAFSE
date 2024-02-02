%% Demonstration for trajectories of Q, P and their derivatives.
% Here we fix the initial condition and potential functionl, i.e., 
%    H(q, p) = (|p|^2+delta^2)^(alpha/2) / alpha + V(q)
%    V(q) = |q|^2 /2

%% Compute derivative in the righthand sides of ODEs 
% syms Q P delta alpha
% T = (P^2+delta^2)^(alpha/2) / alpha;
% V = Q^2/2;
% D1T = diff(T,'P');
% D2T = diff(D1T,'P');
% D3T = diff(D2T,'P');

folder = './figures';
if ~exist(folder, 'file')
    mkdir(folder);
end

%% Parameters setting
global alpha delta
alpha = 1.5;
delta = 0.01*2^-3;

dt = 1e-3;
Tf = 8;
Nt = floor(Tf/dt + 0.1);

t  = zeros(1, Nt+1);
t(1) = 0;

%% Compute the trajectory Q & P
Q0 = 1;
P0 = 1.0;
trj = zeros(2, Nt+1);
trj(:, 1)=[Q0; P0];
k = 1;
while t(k) <= Tf
    dt = min((P0.^2+delta.^2), 1e-2);
    Q     = Q0;
    P     = P0;
    [FQ0,FP0] = FunRight(Q,P);
    %!# 1-step
    Q     = Q0   + dt/2.0*FQ0;
    P     = P0   + dt/2.0*FP0;
    [FQ1,FP1] = FunRight(Q,P);
    %!# 2-step
    Q     = Q0   + dt/2.0*FQ1;
    P     = P0   + dt/2.0*FP1;
    [FQ2,FP2] = FunRight(Q,P);
    %!# 3-step
    Q     = Q0   + dt*FQ2;
    P     = P0   + dt*FP2;
    [FQ3,FP3] = FunRight(Q,P);
    %!# evolution
    Q0    = Q0    + dt * ( FQ0   + 2.0*FQ1   + 2.0*FQ2   + FQ3   ) /6.0;
    P0    = P0    + dt * ( FP0   + 2.0*FP1   + 2.0*FP2   + FP3   ) /6.0;
    %
    trj(:, k+1) = [Q0; P0];
    t(k+1) = t(k) + dt;
    k = k + 1;
end
Nt = k-1;
t = t(1: Nt+1);
dtset = diff(t);
trj = trj(:, 1:Nt+1);
Qfun = spline(t, trj(1, :));
Pfun = spline(t, trj(2, :));
tt = linspace(0, Tf, 30000);
QQ = ppval(Qfun, tt);
PP = ppval(Pfun, tt);

figure;
% plot(t, trj, tt, [QQ; PP]);  % verification
plot(tt, [QQ; PP]);
legend('Q', 'P');
saveas(gcf, './figures/QPfun.png');


%% compute DzQ & DzP
D1Q0 = 1;
D1P0 = -1i;
trj1 = zeros(2, Nt+1);
trj1(:, 1) = [D1Q0; D1P0];
for k = 1 : Nt
    dt = dtset(k);
    %!# 0-step
    t0 = t(k);
    Q = ppval(Qfun, t0);
    P = ppval(Pfun, t0);
    D1Q     = D1Q0;
    D1P     = D1P0;
    [FQ0, FP0] = F1(Q, P, D1Q, D1P);
    %!# 1-step
    t0 = t(k) + dt/2;
    Q = ppval(Qfun, t0);
    P = ppval(Pfun, t0);
    D1Q     = D1Q0   + dt/2.0*FQ0;
    D1P     = D1P0   + dt/2.0*FP0;
    [FQ1, FP1] = F1(Q, P, D1Q, D1P);
    %!# 2-step
    D1Q     = D1Q0   + dt/2.0*FQ1;
    D1P     = D1P0   + dt/2.0*FP1;
    [FQ2, FP2] = F1(Q, P, D1Q, D1P);
    %!# 3-step
    t0 = t(k) + dt;
    Q = ppval(Qfun, t0);
    P = ppval(Pfun, t0);
    D1Q     = D1Q0   + dt*FQ2;
    D1P     = D1P0   + dt*FP2;
    [FQ3, FP3] = F1(Q, P, D1Q, D1P);
    %!# evolution
    D1Q0    = D1Q0    + dt * ( FQ0   + 2.0*FQ1   + 2.0*FQ2   + FQ3   ) /6.0;
    D1P0    = D1P0    + dt * ( FP0   + 2.0*FP1   + 2.0*FP2   + FP3   ) /6.0;
    %
    trj1(:, k+1) = [D1Q0; D1P0];
end
DQfun = spline(t, trj1(1,:));
DPfun = spline(t, trj1(2,:));
D1QQ = ppval(DQfun, tt);
D1PP = ppval(DPfun, tt);

figure;
% plot(t, trj1, tt, [D1QQ; D1PP]);  % verification
plot(tt, real([D1QQ; D1PP]));
legend('DQ', 'DP');
saveas(gcf, './figures/DQDP_real.png')
figure;
plot(tt, imag([D1QQ; D1PP]));
legend('DQ', 'DP');
saveas(gcf, './figures/DQDP_imag.png')


%% compute D2Q & D2P
D2Q0 = 0;
D2P0 = 0;
trj2 = zeros(2, Nt+1);
trj2(:, 1) = [D2Q0; D2P0];
for k = 1 : Nt
    dt = dtset(k);
    %!# 0-step
    t0 = t(k);
    % Q = ppval(Qfun, t0);
    P = ppval(Pfun, t0);
    % D1Q = ppval(DQfun, t0);
    D1P = ppval(DPfun, t0);
    D2Q     = D2Q0;
    D2P     = D2P0;
    [FQ0, FP0] = F2(Q, P, D1Q, D1P, D2Q, D2P);
    %!# 1-step
    t0 = t(k) + dt/2;
    % Q = ppval(Qfun, t0);
    P = ppval(Pfun, t0);
    % D1Q = ppval(DQfun, t0);
    D1P = ppval(DPfun, t0);
    D2Q     = D2Q0   + dt/2.0*FQ0;
    D2P     = D2P0   + dt/2.0*FP0;
    [FQ1, FP1] = F2(Q, P, D1Q, D1P, D2Q, D2P);
    %!# 2-step
    D2Q     = D2Q0   + dt/2.0*FQ1;
    D2P     = D2P0   + dt/2.0*FP1;
    [FQ2, FP2] = F2(Q, P, D1Q, D1P, D2Q, D2P);
    %!# 3-step
    t0 = t(k) + dt;
    % Q = ppval(Qfun, t0);
    P = ppval(Pfun, t0);
    % D1Q = ppval(DQfun, t0);
    D1P = ppval(DPfun, t0);
    D2Q     = D2Q0   + dt*FQ2;
    D2P     = D2P0   + dt*FP2;
    [FQ3, FP3] = F2(Q, P, D1Q, D1P, D2Q, D2P);
    %!# evolution
    D2Q0    = D2Q0    + dt * ( FQ0   + 2.0*FQ1   + 2.0*FQ2   + FQ3   ) /6.0;
    D2P0    = D2P0    + dt * ( FP0   + 2.0*FP1   + 2.0*FP2   + FP3   ) /6.0;
    %
    trj2(:, k+1) = [D2Q0; D2P0];
end

figure;
plot(t, abs(trj2));
legend('D2Q', 'D2P');
saveas(gcf, './figures/D2QD2P_abs.png')
figure;
plot(t, real(trj2));
legend('D2Q', 'D2P');
saveas(gcf, './figures/D2QD2P_real.png')
figure;
plot(t, imag(trj2));
legend('D2Q', 'D2P');
saveas(gcf, './figures/D2QD2P_imag.png')


% ----------------------------------------------------------------
% Subfunctions for the righthand side of ODEs
% ----------------------------------------------------------------

function [FQ, FP] = FunRight(Q, P)
% function to define the FGA ODEs
    global alpha delta
    DT =  P .* (P.^2+delta^2).^(alpha/2-1);
    DV = Q;
    FQ = DT;
    FP = -DV;
end

function [FD1Q, FD1P] = F1(Q, P, D1Q, D1P)
    global alpha delta
    Pdelta = sqrt(P^2+delta^2);
    D2T = Pdelta^(alpha-2) + (alpha-2)*P^2*Pdelta^(alpha-4);
    D2V = 1;
    FD1Q = D1P*D2T;
    FD1P = -D1Q * D2V;
end

function [FD2Q, FD2P] = F2(Q, P, D1Q, D1P, D2Q, D2P)
    global alpha delta
    Pdelta = sqrt(P^2+delta^2);
    D2T = Pdelta^(alpha-2) + (alpha-2)*P^2*Pdelta^(alpha-4);
    D3T = (alpha-2)*Pdelta^(alpha-4)*P + 2*(alpha-2)*Pdelta^(alpha-4)*P+(alpha-2)*Pdelta^(alpha-6)*P^3;
    D2V = 1;
    % D3V=0;
    FD2Q = D2P*D2T + D1P^2 * D3T;
    FD2P = - D2Q * D2V;
    %FD2P  = - D2Q * D2V - D1Q * D3V * D1Q';
end
