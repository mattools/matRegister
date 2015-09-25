function installOolip(varargin)
%INSTALLOOLIP Add the paths required for using OOLIP library
%
%   usage: 
%       installOolip;
%   All the required directories are successively added, in appropriate
%   order to comply with dependencies.
%
%   Example
%   installOolip
%
%   See also
%

% ------
% Author: David Legland
% e-mail: david.legland@grignon.inra.fr
% Created: 2008-01-17,    using Matlab 7.4.0.287 (R2007a)
% Copyright 2008 INRA - BIA PV Nantes - MIAJ Jouy-en-Josas.
% Licensed under the terms of the modified BSD licence

% extract library path
fileName = mfilename('fullpath');
libDir = fileparts(fileName);

moduleNames = {...
    'optimizers', ...
    'transforms', ...
    'interpolators', ...
    'imageSamplers', ...
    'transforms', ...
    'utils', };

disp('Installing OOLIP Library');
addpath(libDir);

% add all library modules
for i = 1:length(moduleNames)
    name = moduleNames{i};
    fprintf('Adding module: %-20s', name);
    addpath(fullfile(libDir, name));
    disp(' (done)');
end

