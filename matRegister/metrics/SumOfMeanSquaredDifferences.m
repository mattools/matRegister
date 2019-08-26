classdef SumOfMeanSquaredDifferences < CostFunction
% Compute the sum of MSD on each image pair.
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
    
end % helper methods


methods
    function varargout = evaluate(obj, params)
        % Evaluate the value and eventually the gradient of the function.
        %
        
        % number of images
        nImg = length(obj.Images);
        
        % ensure we have initialized transformed images and gradients
        if length(obj.TransformedImages)~=nImg
            createTransformedImages(obj);
        end
        if length(obj.TransformedGradients)~=nImg
            createTransformedGradients(obj);
        end
        
        % update transform parameters
        ind1 = 1;
        for i = 1:nImg
            % number of parameters of current transform
            transfo = obj.Transforms{i};
            nParams = getParameterLength(transfo);
            
            % extract parameters corresponding to current transform
            ind2 = ind1 + nParams-1;
            transfoParams = params(ind1:ind2);
            setParameters(transfo, transfoParams);
            
            % update for next transform
            ind1 = ind2 + 1;
        end
        
        % compute metric value and eventually gradient
        if nargout <= 1
            varargout = {computeValue(obj)};
        else
            [fval, grad] = computeValueAndGradient(obj);
            varargout = {fval, grad};
        end
    end
    
    function fval = computeValue(obj)
        % Compute metric value using current state.
        %
        
        % extract number of images
        nImages = length(obj.Images);
        
        % generate all possible couples of images
        combis  = sortrows(combnk(1:nImages, 2));
        nCombis = size(combis, 1);
        
        % compute SSD for each image couple
        res = zeros(nCombis, 1);
        for i = 1:nCombis
            i1 = combis(i,1);
            i2 = combis(i,2);
            
            res(i) = computeMeanSquaredDifferences(...
                obj.Images{i1}, obj.Images{i2}, obj.Points);
        end
        
        % sum of SSD computed over couples
        fval = sum(res);
        
        
        function res = computeMeanSquaredDifferences(img1, img2, points)
            
            % compute values in image 1
            values1 = evaluate(img1, points);
            
            % compute values in image 2
            values2 = evaluate(img2, points);
            
            
            % compute squared differences
            diff = (values2 - values1).^2;
            
            % Sum of squared differences normalized by number of test points
            res = mean(diff);
        end
    end

    
    function [fval, grad] = computeValueAndGradient(obj)
        % Compute metric value and gradient using current state.
        %
        %   [FVAL, GRAD] = computeValueAndGradient(METRIC)
      
        % extract number of images (which is also the number of transforms)
        nImages = length(obj.Images);
        
        % number of parameter for each transform (assumes obj is the same)
        nParams = length(obj.Transforms{1}.Params);
        
        % initialize empty metric value
        fval = 0;
        
        % initialize empty gradient
        grad = zeros(1, nImages*nParams);
        
        % index of current params vactor within global param vector
        ind = 1;
        
        nPoints = size(obj.Points, 1);
        allValues = zeros(nPoints, nImages);
        
        % compute values in each image
        for i = 1:nImages
            [values, inside] = evaluate(obj.TransformedImages{i}, obj.Points);
            allValues(:, i) = values;
            allValues(~inside, i) = 0;
        end
        
        
        % Main iteration over all images (image1)
        % image1 is assumed to be the moving image
        for i = 1:nImages
            % extract data specific to current image
            transfo = obj.Transforms{i};
            gradient = obj.TransformedGradients{i};
            
            % initialize zero gradient vector for current transform
            grad_i = zeros(1, nParams);
            
            % evaluate gradient, and re-compute points within image frame, as
            % gradient evaluator can have different behaviour at image borders.
            [gradVals, gradInside] = evaluate(gradient, obj.Points);
            
            % convert to indices
            inds    = find(gradInside);
            nbInds  = length(inds);
            g = zeros(nbInds, nParams);
            
            for k = 1:nbInds
                iInd = inds(k);
                
                % compute spatial jacobien for current point
                % (in physical coords)
                p0  = obj.Points(iInd, :);
                jac = parametricJacobian(transfo, p0);
                
                % local contribution to metric gradient
                g(iInd, :) = gradVals(iInd, :)*jac;
            end
            
            % second iteration over all other images
            % image2 is fixed image, and we look for average transform towards all
            % images
            for j = [1:i-1 i+1:nImages]
                % compute differences
                diff = allValues(:,i) - allValues(:,j);
                
                % Sum of squared differences normalized by number of test points
                fij = mean(diff.^2);
                
                % re-compute differences, by considering position that can be used
                % for computing gradient
                diff = allValues(gradInside,i) - allValues(gradInside,j);
                
                % compute gradient vectors weighted by local differences
                gd = g(inds,:).*diff(:, ones(1, nParams));
                
                % mean of valid gradient vectors
                gij = mean(gd, 1);
                
                % compute sum of values for all images
                fval    = fval + fij;
                grad_i  = grad_i + gij;
            end
            
            % update global parameter vector
            grad(ind:ind+nParams-1) = grad_i;
            ind = ind + nParams;
        end        
    end
    
end


end % classdef
