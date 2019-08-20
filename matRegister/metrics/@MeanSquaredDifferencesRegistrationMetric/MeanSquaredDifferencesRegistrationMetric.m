classdef MeanSquaredDifferencesRegistrationMetric < RegisteredImageMetric
%MEANSQUAREDDIFFERENCESREGISTRATIONMETRIC  One-line description here, please.
%
%   output = MeanSquaredDifferencesRegistrationMetric(input)
%
%   Example
%   MeanSquaredDifferencesRegistrationMetric
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
    % The gradient image
    GradientImage;
end
 
%% Constructor
methods
    function obj = MeanSquaredDifferencesRegistrationMetric(varargin)
        
        grad = [];
        if nargin==5
            inputs = varargin([1:3 5]);
            grad = varargin{4};
        else
            inputs = varargin;
        end
        
        obj = obj@RegisteredImageMetric(inputs{:});
        obj.GradientImage = grad;
        
    end % constructor
 
end % construction function

%% Some general usage methods
methods
    function transform = getTransform(obj)
        transform = obj.Transform;
    end
    
    function setTransform(obj, transform)
        obj.Transform = transform;
        % build the backward transformed image
        obj.TransformedImage = BackwardTransformedImage(...
            obj.MovingImage, obj.Transform);
    end
    
    function gradient = getGradientImage(obj)
        gradient = obj.GradientImage;
    end
    
    function setGradientImage(obj, gradient)
        obj.GradientImage = gradient ;
    end 
        
end

end % classdef
