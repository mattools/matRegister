classdef ImageToImageMetric < BaseFunction
%IMAGETOIMAGEMETRIC Abstract class that define a metric between 2 images.
%
%   output = ImageToImageMetric(input)
%
%   IMG1 and IMG2 should be instance of ImageFunction. If they are not,
%   they are converted using LinearInterpolator by default.
%
%   Example
%   M = ImageToImageMetric(IMG1, IMG2, PTS)
%   M.evaluate()
%
%   See also
%     MeanSquaredDifferencesMetric
%

% ------
% Author: David Legland
% e-mail: david.legland@inra.fr
% Created: 2010-08-12,    using Matlab 7.9.0.529 (R2009b)
% Copyright 2010 INRA - Cepia Software Platform.

%% Properties
properties
    % the fixed image
    Img1;
    
    % the moving image
    Img2;
    
    % the set of points that will be used for evaluating the metric
    Points;
    
end % properties

%% Constructor
methods (Access = protected)
    
    function obj = ImageToImageMetric(varargin)
        % Protected constructor of ImageToImageMetric class
        errorID = 'ImageToImageMetric:Constructor';
        
        if nargin==2
            % Inputs are IMG1 and IMG2
            error(errorID, 'not yet implemented');
            
        elseif nargin==3
            % Inputs are IMG1, IMG2, and the set of points
            
            % setup image 1
            var = varargin{1};
            if isa(var, 'ImageFunction')
                obj.Img1 = var;
            elseif isa(var, 'Image')
                obj.Img1 = ImageInterpolator.create(var, 'linear');
            elseif isnumeric(var)
                img = Image.create(var);
                obj.Img1 = ImageInterpolator.create(img, 'linear');
            else
                error(errorID, 'First input should be an ImageFunction');
            end
            
            % setup image 2
            var = varargin{2};
            if isa(var, 'ImageFunction')
                obj.Img2 = var;
            elseif isa(var, 'Image')
                obj.Img2 = ImageInterpolator.create(var, 'linear');
            elseif isnumeric(var)
                img = Image.create(var);
                obj.Img2 = ImageInterpolator.create(img, 'linear');
            else
                error(errorID, 'Second input should be an ImageFunction');
            end
            
            % setup points
            obj.points = varargin{3};
            
        else
            error(errorID, 'Need at least 2 inputs');
        end

    end % constructor
    
end % methods

%% Abstract methods
methods (Abstract)
    computeValue(obj)
    % Evaluate the difference between the two images.
end

end % classdef
