function FGAFSE_startup()
% FGAFSE_STARTUP Startup file for FGAFSE
%   Add paths of FGAFSE to Matlab.

%  Copyright (c) 2024 Hengzhun Chen, Fudan University,
%                     Lihui Chai, Sun Yat-sen University.
%  This file is distributed under the terms of the MIT License.


file_path = mfilename('fullpath');
tmp = strfind(file_path, 'FGAFSE_startup');
file_path = file_path(1 : (tmp(end)-1));

% Folder for all source files recursively
addpath(genpath([file_path 'src']));

% Folder for all test files recursively
addpath(genpath([file_path 'test']));

end