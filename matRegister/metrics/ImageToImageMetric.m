classdef ImageToImageMetric < BaseFunction
% Abstract class that define a metric between 2 images.
%
%   M = ImageToImageMetric(IMG1, IMG2, PTS)
%   IMG1 and IMG2 should be instance of ImageFunction. If they are not,
%   they are converted using LinearInterpolator by default.
%
%   Example
%   M = ImageToImageMetric(IMG1, IMG2, PTS)
%   evaluate(M)
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

methods
    function setFixedImage(this, img)
        % Setup fixed image of registration algorithm.
        %
        %   setFixedImage(METRIC, IMG)
        %   IMG can be either an instance of image, or an instance of
        %   ImageFunction, for example an interpolated image.
        
        this.Img1 = img;
    end
    
    function setMovingImage(this, img)
        % Setup moving image of registration algorithm.
        %
        %   setMovingImage(METRIC, IMG)
        %   IMG can be either an instance of image, or an instance of
        %   ImageFunction, for example an interpolated image.
        
        this.Img2 = img;
    end
    
    function setPoints(this, points)
        % Setup the inner set of test points for this metric.
        %
        %   setPoints(METRIC, POINTS);
        this.Points = points;
    end
end

end % classdef
