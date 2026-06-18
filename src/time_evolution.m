function [A, S, Q, P, DzQ, DzP] = time_evolution( ...
        A0, S0, Q0, P0, DzQ0, DzP0, dt, finalTime, odes, solver, evalKinetic, evalPotential)
% TIME_EVOLUTION Evolve the FGA ODEs with the selected time integrator.
%    Inputs:
%        A0, S0, Q0, P0, DzQ0, DzP0
%                  -- initial FGA parameter arrays of size nGB
%        dt        -- maximum step size of time t
%        finalTime -- duration of time evolution
%        odes      -- closed function handle of the FGA ODEs
%                     [DQ, DP, DS, DlogA, DtDzQ, DtDzP] = ...
%                             odes(Q, P, DzQ, DzP)
%        solver    -- ODE solver: 'rk4' (default) or 'symplectic'
%        evalKinetic
%                  -- required for 'symplectic'; closed function handle
%                     [T, DT, D2T] = evalKinetic(P)
%        evalPotential
%                  -- required for 'symplectic'; closed function handle
%                     [V, DV, D2V] = evalPotential(Q)
%                     It should adapt the dimension-specific user potential
%                     to accept Q as an nGB-by-d array.
%    Outputs:
%        A, S, Q, P, DzQ, DzP
%                  -- FGA parameter arrays of size nGB
%    DzQ and DzP are returned for diagnostics or continued evolution.
%
%    See also TIME_EVOLUTION_RK4, TIME_EVOLUTION_SYMPLECTIC, FGA1d, FGA2d.

%  Copyright (c) 2024 Hengzhun Chen, Fudan University,
%                     Lihui Chai, Sun Yat-sen University.
%  This file is distributed under the terms of the MIT License.


if nargin < 10 || isempty(solver)
    solver = 'rk4';
end

solver = validatestring(solver, {'rk4', 'symplectic', 'stormer-verlet'});
if strcmp(solver, 'stormer-verlet')
    solver = 'symplectic';
end

switch solver
    case 'rk4'
        [A, S, Q, P, DzQ, DzP] = time_evolution_rk4( ...
            A0, S0, Q0, P0, DzQ0, DzP0, dt, finalTime, odes);
    case 'symplectic'
        if nargin < 12 || isempty(evalKinetic) || isempty(evalPotential)
            error('FGAFSE:MissingHamiltonianData', ...
                'The symplectic solver requires kinetic and potential data functions.');
        end
        [A, S, Q, P, DzQ, DzP] = time_evolution_symplectic( ...
            A0, S0, Q0, P0, DzQ0, DzP0, dt, finalTime, evalKinetic, evalPotential);
end

end
