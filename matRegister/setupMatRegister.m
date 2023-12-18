function setupMatRegister(varargin)
% Add the paths required for using MatRegister library.
%
%   usage: 
%       setupMatRegister;
%   All the required directories are successively added, in the appropriate
%   order to comply with dependencies.
%
%   Example
%   setupMatRegister
%
%   See also
%

% ------
% Author: David Legland
% e-mail: david.legland@inrae.fr
% Created: 2008-01-17,    using Matlab 7.4.0.287 (R2007a)
% Copyright 2008 INRAE - BIA PV Nantes - MIAJ Jouy-en-Josas.
% Licensed under the terms of the modified BSD licence

% extract library path
fileName = mfilename('fullpath');
libDir = fileparts(fileName);

moduleNames = {...
    'imageSamplers', ...
    'transforms', ...
    'interpolators', ...
    'metrics', ...
    'shapeMetrics', ...
    'optimizers', ...
    'utils', };

disp('Installing MatRegister Library');
addpath(libDir);

% add all library modules
for i = 1:length(moduleNames)
    name = moduleNames{i};
    fprintf('Adding module: %-20s', name);
    addpath(fullfile(libDir, name));
    disp(' (done)');
end

