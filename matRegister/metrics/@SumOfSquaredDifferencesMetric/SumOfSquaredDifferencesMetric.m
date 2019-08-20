classdef SumOfSquaredDifferencesMetric < ImageToImageMetric
%SumOfSquaredDifferencesMetric Compute sum of squared differences metric
%
%   output = SumOfSquaredDifferencesMetric(input)
%
%   Example
%   SumOfSquaredDifferencesMetric
%
%   See also
%

% ------
% Author: David Legland
% e-mail: david.legland@inra.fr
% Created: 2010-08-12,    using Matlab 7.9.0.529 (R2009b)
% Copyright 2010 INRA - Cepia Software Platform.

%% Some properties specific to the SSD metric
properties
    % the transform model (necessary to compute the metric gradient)
    Transform;
    
    % The gradient image of the moving image (necessary to compute the metric gradient)
    GradientImage;
end

%% Constructor
methods
    function obj = SumOfSquaredDifferencesMetric(varargin)
        % calls the parent constructor
        obj = obj@ImageToImageMetric(varargin{:});
        
    end % constructor
    
end % methods

%% Accessors and modifiers
methods
    function transform = getTransform(obj)
        transform = obj.Transform;
    end
    
    function setTransform(obj, transform)
        obj.Transform = transform;
    end
    
    function gradient = getGradientImage(obj)
        gradient = obj.GradientImage;
    end
    
    function setGradientImage(obj, gradient)
        obj.GradientImage = gradient ;
    end 
        
end % methods

end % classdef
