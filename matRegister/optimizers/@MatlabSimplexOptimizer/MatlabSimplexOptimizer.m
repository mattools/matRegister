classdef MatlabSimplexOptimizer < Optimizer
%MATLABSIMPLEXOPTIMIZER Encapsulation of Matlab function fminsearch
%
%   output = MatlabSimplexOptimizer(input)
%
%   Example
%   % Run the simplex otpimizer on the Rosenbrock function
%     optimizer = MatlabSimplexOptimizer;
%     optimizer.setCostFunction(@rosenbrock);
%     optimizer.setParameters([0 0]);
%     [xOpt value] = optimizer.startOptimization();
%     xOpt
%       xOpt =
%           0.9993    0.9984
%     value
%       value =
%           5.5924e-006
%
%
%   Requires
%   Optimization Toolbox
%
%   See also
%     NelderMeadSimplexOptimizer
%

% ------
% Author: David Legland
% e-mail: david.legland@inra.fr
% Created: 2010-10-06,    using Matlab 7.9.0.529 (R2009b)
% Copyright 2010 INRA - Cepia Software Platform.

%% Properties
properties
    NIters = 200;
end

%% Constructor
methods
    function obj = MatlabSimplexOptimizer(varargin)
        obj = obj@Optimizer();
    end
end

end