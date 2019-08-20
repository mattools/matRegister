classdef ParameterizedFunctionEvaluator < CostFunction
% Update parameters of a transform and evaluate it.
%
%   output = ParameterizedFunctionEvaluator(input)
%
%   Example
%   ParameterizedFunctionEvaluator
%
%   See also
%

% ------
% Author: David Legland
% e-mail: david.legland@inra.fr
% Created: 2010-11-03,    using Matlab 7.9.0.529 (R2009b)
% Copyright 2010 INRA - Cepia Software Platform.


%% Properties
properties
    Transform;
    BaseFunction;
end


%% Constructor
methods
    function obj = ParameterizedFunctionEvaluator(transform, baseFunction)
        %Constructor for a new MetricEvaluator
        % THIS = MetricEvaluator(TRANSFO, METRIC)
        
        % test class of transform
        if ~isa(transform, 'ParametricObject')
            error('First argument must be a parametric object');
        end
        
        obj.Transform = transform;
        obj.BaseFunction = baseFunction;
      
    end % constructor
 
end % construction function
 

%% General methods
methods
 
    function [res, grad] = evaluate(obj, params)
        % Update transform parameters, and compute function value (and grad)
        
        % update params
        setParameters(obj.Transform, params);
        
        % compute value and eventually graident
        if nargout <= 1
            res = computeValue(obj.BaseFunction);
        else
            [res, grad] = computeValueAndGradient(obj.BaseFunction);
        end
    end
end % general methods
 
end % classdef

