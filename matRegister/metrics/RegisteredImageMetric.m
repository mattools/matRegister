classdef RegisteredImageMetric < handle
%REGISTEREDIMAGEMETRIC  One-line description here, please.
%
%   output = RegisteredImageMetric(input)
%
%   Example
%   RegisteredImageMetric
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
    % the reference image (interpolated)
    FixedImage;

    % the moving image (interpolated)
    MovingImage;

    % the transform model from fixed image to moving image
    Transform;

    % the set of points used for metric evaluation
    Points;

    % the transformed image
    TransformedImage;
end

%% Constructor
methods (Access = protected)
    function obj = RegisteredImageMetric(img1, img2, transfo, points)
        
        % check inputs
        if nargin<3
            error('Need to specify some input arguments...');
        end
        
        % process fixed image
        if isa(img1, 'ImageFunction')
            obj.FixedImage = img1;
        else
            obj.FixedImage = ImageInterpolator.create(img1, 'linear');
        end
        
        % process moving image
        if isa(img2, 'ImageFunction')
            obj.MovingImage = img2;
        else
            obj.MovingImage = ImageInterpolator.create(img2, 'linear');
        end
        
        % in the case of 3 arguments, the third one is the number of points
        if nargin == 3
            obj.Points = transfo;
            % create an identify transform
            nd = ndims(img1);
            obj.Transform = TranslationModel(nd);
            
        else
            % process transform
            if ~isa(transfo, 'ParametricTransform')
                error('Transform should be parametric');
            end
            obj.Transform = transfo;
            
            % stores the set of test points
            obj.Points = points;
        end
        
        % build the backward transformed image
        obj.TransformedImage = BackwardTransformedImage(...
            obj.MovingImage, obj.Transform);
        
    end % constructor

end % construction function

%% Abstract methods
methods (Abstract)
    res = evaluate(obj, params)
    % Evaluate the image to image metric
    
end % general methods

end % classdef
