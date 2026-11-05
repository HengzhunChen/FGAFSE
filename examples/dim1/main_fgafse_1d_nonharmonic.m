%% 1D FGA error study with a non-harmonic double-well potential
%
% This follows main_fgafse_1d.m, but replaces the harmonic trap by
% V(x) = ((x-4)^2-1)^2/16.  The two delta scalings requested in
% docs/new_numerical_concern.md are run and saved separately.

exampleDir = fileparts(mfilename('fullpath'));
cd(fullfile(exampleDir, '..', '..'));
FGAFSE_startup();
cd(exampleDir);

alpha = [1.1, 1.3, 1.5, 1.7, 1.9];
alphaRegularized = [1.75, 1.8, 1.85, 1.9, 1.95];
vepsExp = [-6, -7, -8, -9, -10];
right_x = 8;
final_time = 5;
useCachedTSSA = false; % The potential differs from the harmonic-example cache.
fgaDt = 5e-4; % Change this value to study FGA ODE time-step sensitivity.

studies = struct('power', {1, 7/12}, ...
    'tag', {'delta_eps', 'delta_eps_7div12'}, ...
    'alpha', {alpha, alphaRegularized}, ...
    'figure', {'L2err_nonharmonic_delta_eps.eps', ...
               'L2err_nonharmonic_delta_eps_7div12.eps'});

results = struct();
for studyIndex = 1:numel(studies)
    study = studies(studyIndex);
    alphaStudy = study.alpha;
    fprintf('\nNon-harmonic study: delta = epsilon^%g\n', study.power);
    figureName = fullfile(exampleDir, study.figure);
    error_decay_1d(alphaStudy, vepsExp, right_x, final_time, @initWave, ...
        @potential, figureName, useCachedTSSA, study.power, fgaDt);

    % error_decay_1d writes its numerical table to L2_error.mat. Preserve a
    % separate copy for each delta scaling and add the study metadata.
    saved = load('L2_error.mat');
    saved.alpha = alphaStudy;
    saved.vepsExp = vepsExp;
    saved.right_x = right_x;
    saved.final_time = final_time;
    saved.deltaPower = study.power;
    saved.fgaDt = fgaDt;
    saved.potentialName = 'V(x)=((x-4)^2-1)^2/16';
    results.(study.tag) = saved;
    save(fullfile(exampleDir, ['L2_error_nonharmonic_' study.tag '.mat']), '-struct', 'saved');
end
save(fullfile(exampleDir, 'L2_error_nonharmonic.mat'), '-struct', 'results');

function u0 = initWave(x, veps)
    q0 = 4;
    p0 = 1;
    u0 = (128 / pi)^(1/4) * exp(-64 * (x - q0).^2) ...
        .* exp(1i * p0 * x / veps);
end

function [V, DV, D2V] = potential(Q)
    y = Q - 4;
    V = (y.^2 - 1).^2 / 16;
    DV = y .* (y.^2 - 1) / 4;
    D2V = (3 * y.^2 - 1) / 4;
end
