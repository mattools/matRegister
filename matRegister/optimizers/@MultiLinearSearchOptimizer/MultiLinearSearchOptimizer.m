classdef MultiLinearSearchOptimizer < Optimizer
%MultiLinearSearchOptimizer Optimize along successive directions
%
%   output = MultiLinearSearchOptimizer(input)
%
%   Example
%   % Run the simplex otpimizer on the Rosenbrock function
%     optimizer = MultiLinearSearchOptimizer;
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
%   See also
%

% ------
% Author: David Legland
% e-mail: david.legland@inra.fr
% Created: 2010-10-06,    using Matlab 7.9.0.529 (R2009b)
% Copyright 2010 INRA - Cepia Software Platform.

%% Properties
properties
    NIters = 200;
    
    DirectionSet;
    
    % TODO: change to another name/type of parameter
    Verbose =true;
end

%% Constructor
methods
    function obj = MultiLinearSearchOptimizer(varargin)
        obj = obj@Optimizer(varargin{:});
    end
end

%% General methods
methods
    function directionSet = getDirectionSet(obj)
        directionSet = obj.DirectionSet;
    end
    
    function setDirectionSet(obj, newDirectionSet)
        obj.DirectionSet = newDirectionSet;
    end
end

%% Private methods
methods (Access = private)
    function initDirectionSet(obj)
        % init a new set of direction to span the paramater space
        
        % allocate memory
        nParams = length(obj.Params);
        set = zeros(nParams+ceil(nParams/2), nParams);
        
        % init main directions (isothetic directions)
        set(1:nParams, 1:nParams) = eye(nParams);
        
        % init also diagonal directions
        tmp = eye(nParams) + diag(ones(nParams-1, 1), 1);
        tmp = tmp(1:2:end, :);
        set(nParams+1:end, :) = tmp;

        % process last line
        if floor(nParams/2) ~= ceil(nParams/2)
            set(end, 1) = 1;
        end
        
        % stores result
        obj.DirectionSet = set;
    end
end

end