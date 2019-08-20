classdef GaussianLinearSearchOptimizer < Optimizer
%GaussianLinearSearchOptimizer Optimize each parameter individually
%
%   output = GaussianLinearSearchOptimizer(input)
%
%   Example
%     x0 = [0 0];
%     optim = GaussianLinearSearchOptimizer(@rosenbrock, x0);
%     optim.nValues = 50;
%     optim.nIter = 100;
%     optim.parameterVariability = .2;
%     xHat = optim.startOptimization([0 0])
%     xHat =
%         0.5469    0.2951
%
%   See also
%

% ------
% Author: David Legland
% e-mail: david.legland@inra.fr
% Created: 2010-11-24,    using Matlab 7.9.0.529 (R2009b)
% Copyright 2010 INRA - Cepia Software Platform.


properties
    % the 'variability' of each parameter. Default is 1.
    ParameterVariability;
    
    % number of values to compute on each param
    NValues = 50;
    
    % maximum number of iterations
    NIters = 10;
    
end

methods
    function obj = GaussianLinearSearchOptimizer(varargin)
        obj = obj@Optimizer(varargin{:});
    end
end

end