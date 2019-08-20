classdef RegularizedImageToImageMetric < ParametricFunction
%REGULARIZEDIMAGETOIMAGEMETRIC  One-line description here, please.
%
%   output = RegularizedImageToImageMetric(input)
%
%   Example
%   RegularizedImageToImageMetric
%
%   See also
%

% ------
% Author: David Legland
% e-mail: david.legland@inra.fr
% Created: 2010-10-27,    using Matlab 7.9.0.529 (R2009b)
% Copyright 2010 INRA - Cepia Software Platform.


 
%% Properties
properties
    % an image to image metric
    Metric;
    
    % the regularization object of the transform
    Regul;
    
    % the transform itself
    Transform;
    
    % the coefficient associated to the metric
    Alpha;
    
    % the coefficient associated to the regularisation
    Beta;
end
 
%% Constructor
methods
    function obj = RegularizedImageToImageMetric(metric, regul, alpha, beta)
        obj.Metric = metric;
        obj.Regul  = regul;
        obj.Alpha  = alpha;
        obj.Beta   = beta;
        
    end % constructor
 
end % construction function
 
%% General methods
methods
    function setParameters(obj, params)
        setParameters(obj.Transform, params);
    end
    
    function res = computeValue(obj)
        res1 = obj.Alpha   * computeValue(obj.Metric);
        res2 = obj.Beta    * computeValue(obj.Regul);
        res  = res1 + res2;
    end
    
end % general methods
 
end % classdef

