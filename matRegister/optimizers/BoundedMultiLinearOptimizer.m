classdef BoundedMultiLinearOptimizer < Optimizer
% Optimize each parameter individually.
%
%   output = BoundedMultiLinearOptimizer(input)
%
%   Example
%   BoundedMultiLinearOptimizer
%
%   See also
%

% ------
% Author: David Legland
% e-mail: david.legland@inra.fr
% Created: 2010-11-19,    using Matlab 7.9.0.529 (R2009b)
% Copyright 2010 INRA - Cepia Software Platform.


%% Properties

properties
    % the bounds for each parameter. Should be a Np-by-2 array.
    Bounds;
    
    % number of values to compute on each param
    NValues = 50;
    
    % maximum number of iterations
    NIters = 10;
    
end


%% Constructor

methods
    function obj = BoundedMultiLinearOptimizer(varargin)
        obj = obj@Optimizer(varargin{:});
    end
end


%% Methods from Optimizer interface

methods
    function [params, value] = startOptimization(obj)
        % Run the optimizer, and return optimized parameters.
        %
        %   PARAMS = startOptimization(OPTIM)
        %   PARAMS = OPTIM.startOptimization()
        %
        %   [PARAMS VALUE] = startOptimization(OPTIM)
        %   [PARAMS VALUE] = OPTIM.startOptimization()
        %
        %   Example
        %   StartOptimization
        %
        %   See also
        %
        
        % ------
        % Author: David Legland
        % e-mail: david.legland@inra.fr
        % Created: 2010-11-19,    using Matlab 7.9.0.529 (R2009b)
        % Copyright 2010 INRA - Cepia Software Platform.
        
        
        % Notify beginning of optimization
        notify(obj, 'OptimizationStarted');
        
        params = obj.Params;
        if ~isempty(obj.InitialParameters)
            params = obj.InitialParameters;
        end
        
        nParams = length(params);
        bounds = obj.Bounds;
        
        if size(bounds, 1)~=nParams || size(bounds, 2)~=2
            warning('oolip:NonInitializedParameter',...
                'Bounds not specified, use default bounds');
            bounds = [params'-10 params'+10];
        end
        
        
        for i = 1:obj.NIters
            for p = 1:nParams
                % choose a set of values
                par0 = bounds(p, 1);
                par1 = bounds(p, 2);
                paramValues = linspace(par0, par1, obj.NValues);
                
                % compute the metric for the given param
                res = zeros(obj.NValues, 1);
                for k = 1:obj.NValues
                    params(p) = paramValues(k);
                    res(k) = obj.CostFunction(params);
                end
                
                [value, bestK] = min(res);
                params(p) = paramValues(bestK);
                
                % update optimizer internal state
                obj.Params = params;
                obj.Value  = value;
                
                % notify
                notify(obj, 'OptimizationIterated');
            end
        end
        
        % update inner data
        obj.Params = params;
        obj.Value = value;
        
        % Notify the end of optimization
        notify(obj, 'OptimizationTerminated');
    end
end

end