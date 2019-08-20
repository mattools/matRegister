classdef SumOfMeanSquaredDifferences < CostFunction
%SUMOFMEANSQUAREDDIFFERENCES Compute the sum of MSD on each image couples
%
%   Quite similar to class "SumOfSSDImageSetMetric", but obj class is
%   intended to provide computation of Gradient with respect to parameters.
%
%   output = SumOfMeanSquaredDifferences(input)
%
%   Example
%   SumOfMeanSquaredDifferences
%
%   See also
%

% ------
% Author: David Legland
% e-mail: david.legland@inra.fr
% Created: 2011-01-06,    using Matlab 7.9.0.529 (R2009b)
% Copyright 2011 INRA - Cepia Software Platform.

%% Properties
properties
    % First, some "public" data.
    
    % the set of images to transform
    Images;
    
    % the set of test points in reference space
    Points;
    
    % the set of parametric transforms
    Transforms;
    
    % the set of gradient image. Either as a set of Vector Images, or as a
    % set of Gradient evaluators
    Gradients;
    
    % Then, some private data used during computation
    
    % images after transforms
    TransformedImages;
    
    % gradient images after transforms
    TransformedGradients
    
end 

%% Constructor
methods
    function obj = SumOfMeanSquaredDifferences(images, points, transforms, varargin)
        % Construct a new image set metric
        %
        % Metric = SumOfMeanSquaredDifferences(IMGS, PTS, TRANSFOS);
        % Metric = SumOfMeanSquaredDifferences(IMGS, PTS, TRANSFOS, GRADS);
        %
        
        errorID = 'SumOfMeanSquaredDifferences:Constructor';
        
        if nargin >= 3
            % Inputs are the set of images, the array of points, and the
            % set of transforms
            obj.Images     = images;
            obj.Points     = points;
            obj.Transforms = transforms;
        else
            error(errorID, 'Need at least 3 inputs');
        end
        
        % Gradients can be specified
        if nargin> 3 
            obj.Gradients = varargin{1};
        else
            error('Not yet implemented');
            %TODO: create gradients or gradient interpolators
        end
        
        ensureImagesValidity(obj);
        createTransformedImages(obj);

        ensureGradientsValidity(obj);
        createTransformedGradients(obj);
        
    end % constructor 

end % construction function

%% General methods
methods (Access = protected)

    function ensureImagesValidity(obj)
        % Checks that images are instances of ImageFunction, otherwise
        % creates appropriate interpolators
        for i = 1:length(obj.Images)
            img = obj.Images{i};
            
            % input image should be an image function
            if isa(img, 'ImageFunction')
                continue;
            end
            
            % Try to create an appropriate image function
            if isa(img, 'Image')
                img = ImageInterpolator.create(img, 'linear');
            elseif isnumeric(img)
                img = Image.create(img);
                img = ImageInterpolator.create(img, 'linear');
            else
                error('Unknown image data');
            end
            obj.Images{i} = img;
            
        end
    end
    
    function createTransformedImages(obj)
        % Creates instances of transformed images
        
        % small check
        nImg = length(obj.Images);        
        if length(obj.Transforms) ~= nImg
            error('Number of transforms should equal number of images');
        end
        
        % create array
        obj.TransformedImages = cell(1, nImg);
        
        % iterate over image-transform couple to create transformed images
        for i = 1:nImg
            image = obj.Images{i};
            transfo = obj.Transforms{i};
            tim = BackwardTransformedImage(image, transfo);
            obj.TransformedImages{i} = tim;
        end
    end
    
    function ensureGradientsValidity(obj)
        % Checks that gradient images are instances of ImageFunction, 
        % otherwise creates appropriate interpolators
        for i = 1:length(obj.Gradients)
            gradImg = obj.Gradients{i};
            
            % input image should be an image function
            if isa(gradImg, 'ImageFunction')
                continue;
            end
            
            % Try to create an appropriate image function
            if isa(gradImg, 'Image')
                gradImg = ImageInterpolator.create(gradImg, 'nearest');
            else
                error('Unknown image data');
            end
            obj.Gradients{i} = gradImg;
            
        end
    end
    
    function createTransformedGradients(obj)
        % Creates instances of transformed gradient images or functions
        
        % small check
        nImg = length(obj.Gradients);        
        if length(obj.Transforms) ~= nImg
            error('Number of transforms should equal number of gradient images');
        end
        
        % create array
        obj.TransformedGradients = cell(1, nImg);
        
        % iterate over image-transform couple to create transformed images
        for i = 1:nImg
            image = obj.Gradients{i};
            transfo = obj.Transforms{i};
            tim = BackwardTransformedImage(image, transfo);
            obj.TransformedGradients{i} = tim;
        end
    end
    
end % abstract methods

end % classdef
