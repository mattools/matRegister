classdef MultiLinearSearchOptimizer < Optimizer
% Optimize along successive directions.
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

methods
    function [params, value] = startOptimization(obj)
        %STARTOPTIMIZATION Run the optimization algorithm
        
        
        % ensure there is a valid direction set
        if isempty(obj.DirectionSet)
            initDirectionSet(obj);
        end
        
        
        % default tolerance
        tol = 1e-5;
        
        nDirs = size(obj.DirectionSet, 1);
        
        % main loop
        % Use linear search in a set of directions that span the parameters of the
        % first transform.
        % The same process could be obteained by calling 'directionSetMinimizer'.
        dirIndex = 1;
        for i = 1:obj.NIters
            if obj.Verbose
                disp(sprintf('iteration %d / %d', i, obj.NIters)); %#ok<DSPS>
            end
            
            % current direction
            dir = obj.DirectionSet(dirIndex, :);
            
            % update direction index
            dirIndex = mod(dirIndex, nDirs) + 1;
            
            % use a function handle of 1 variable
            fun1 = @(t) obj.CostFunction(obj.Params + t*dir);
            
            % guess initial bounds of the function
            ax = 0;
            bx = 1;
            [ax, bx, cx] = fMinBracket(fun1, ax, bx);
            
            % search minimum along dimension DIR
            [tmin, value] = brentLineSearch(fun1, ax, bx, cx, tol);
            
            % construct new optimal point
            obj.Params = obj.Params + tmin*dir;
            
            % Call an output function for processing current point
            if ~isempty(obj.OutputFunction)
                % setup optim values
                optimValues.fval = value;
                optimValues.iteration = i;
                optimValues.procedure = 'Linear search';
                
                % call output function with appropriate parameters
                stop = obj.OutputFunction(obj.Params, optimValues, 'iter');
                if stop
                    break;
                end
            end
        end
        
        % format output arguments
        params = obj.Params;
        
        % we need a final end because nested functions are used
    end

end
end