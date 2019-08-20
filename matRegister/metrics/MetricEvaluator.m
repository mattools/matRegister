classdef MetricEvaluator < handle
%METRICEVALUATOR Evaluate a metric that depend on a parametric object
%
%   Usage:
%   EV = MetricEvaluator(TRANSFO, METRIC);
%   EV.evaluate(PARAMS);
%
%   Example
%   MetricEvaluator
%
%   See also
%     ParametricTransform, ImageToImageMetric
%

% ------
% Author: David Legland
% e-mail: david.legland@inra.fr
% Created: 2010-10-29,    using Matlab 7.9.0.529 (R2009b)
% Copyright 2010 INRA - Cepia Software Platform.

 
%% Properties
properties
    Transform;
    Metric;
end


%% Constructor
methods
    function obj = MetricEvaluator(transform, metric)
        %Constructor for a new MetricEvaluator
        % THIS = MetricEvaluator(TRANSFO, METRIC)
        
        % test class of transform
        if ~isa(transform, 'ParametricTransform')
            error('First argument must be a parametric transform');
        end
        
        obj.Transform = transform;
        obj.Metric = metric;
      
    end % constructor
 
end % construction function
 

%% General methods
methods
 
    function [res, grad] = evaluate(obj, params)
        % Update transform parameters, and compute metric value (and grad)
        
        % uupdate params
        setParameters(obj.Transform, params);
        
        % compute value and eventually gradient
        if nargout <= 1
            res = computeValue(obj.Metric);
        else
            [res, grad] = computeValueAndGradient(obj.Metric);
        end
    end
end % general methods
 
end % classdef

