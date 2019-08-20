classdef MeanSquaredDifferencesMetric < ImageToImageMetric
%Image to Image metric that compute mean of squared differences.
%
%   METRIC = MeanSquaredDifferencesMetric(IMG1, IMG2, POINTS)
%   IMG1 and IMG2 are preferentially instances of BackwardTransformedImage,
%   POINTS is a set of ND points stored in NP-by-ND array of double.
%
%   Example
%   MeanSquaredDifferencesMetric
%
%   See also
%

% ------
% Author: David Legland
% e-mail: david.legland@inra.fr
% Created: 2010-08-12,    using Matlab 7.9.0.529 (R2009b)
% Copyright 2010 INRA - Cepia Software Platform.

%% Some properties specific to the MSD metric
properties
    % the transform model (necessary to compute the metric gradient)
    Transform;
    
    % The gradient image
    GradientImage;
end

%% Constructor
methods
    function obj = MeanSquaredDifferencesMetric(varargin)
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
