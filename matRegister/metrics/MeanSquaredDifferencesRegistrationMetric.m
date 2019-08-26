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

methods
    function [res, grad] = evaluate(obj, params)
        % Evaluate value and metric of MSD image to image metric.
        
        setParameters(obj.Transform, params);
        
        if nargout > 1
            [res, grad] = computeValueAndGradient(obj);
        else
            res = computeValue(obj);
        end
    end
    
    function [res, isInside] = computeValue(obj)
        %COMPUTEVALUE Compute metric value.
        %
        % [VALUE, INSIDE] = METRIC.computeValue();
        % Computes and return the value. Returns also a flag that indicates
        % which test points belong to both images.
        %
        
        % compute values in image 1
        [values1, inside1] = evaluate(obj.FixedImage, obj.Points);
        
        % compute values in image 2
        [values2, inside2] = evaluate(obj.TransformedImage, obj.Points);
        
        % keep only valid values
        isInside = inside1 & inside2;
        
        % compute result
        diff = (values2(isInside) - values1(isInside)).^2;
        res = mean(diff);
    end
    
    function [res, grad, isInside] = computeValueAndGradient(obj, varargin)
        % Compute metric value and gradient.
        %
        %   [RES, DERIV] = computeValueAndGradient(METRIC);
        %   This syntax requires that both fields 'transform' and 'gradientImage'
        %   have been initialized.
        %
        %   [RES, DERIV] = computeValueAndGradient(METRIC, TRANSFO, GRADX, GRADY);
        %   [RES, DERIV] = computeValueAndGradient(METRIC, TRANSFO, GRADX, GRADY, GRADZ);
        %   This (deprecated) syntax passes transform model and gradient components
        %   as input arguments.
        %
        % Example:
        % transfo = Translation2DModel([1.2 2.3]);
        % ssdMetric = SumOfSquaredDifferencesMetric(img1, img2, transfo);
        % res = computeValueAndGradient(ssdMetric, model);
        %
        
        if isempty(obj.Transform) || isempty(obj.GradientImage)
            error('Either transform or gradient image was not initialized');
        end
        if isempty(obj.TransformedImage)
            error('Transformed image was not initialized');
        end
        
        nd = ndims(obj.FixedImage);
        if nd == 2
            [res, grad, isInside] = computeValueAndGradientLocal2d(obj);
        elseif nd == 3
            [res, grad, isInside] = computeValueAndGradientLocal3d(obj);
        else
            [res, grad, isInside] = computeValueAndGradientLocal(obj);
        end
        
        
        % end of main function
        
        
        function [res, grad, isInside] = computeValueAndGradientLocal(obj)
            
            % error checking
            if isempty(obj.Transform)
                error('Gradient computation requires transform');
            end
            if isempty(obj.GradientImage)
                error('Gradient computation requires a gradient image');
            end
            
            % compute values in image 1
            [values1, inside1] = evaluate(obj.FixedImage, obj.Points);
            
            % compute values in image 2
            [values2, inside2] = evaluate(obj.TransformedImage, obj.Points);
            
            % keep only valid values
            isInside = inside1 & inside2;
            
            % compute result
            diff = values2(isInside) - values1(isInside);
            res = mean(diff.^2);
            
            %fprintf('Initial SSD: %f\n', res);
            
            % convert to indices
            inds = find(isInside);
            nbInds = length(inds);
            
            transfo = obj.Transform;
            nParams = length(getParameters(transfo));
            g = zeros(nbInds, nParams);
            
            % convert from physical coordinates to index coordinates
            % (assumes spacing is 1 and origin is 0)
            points2 = transformPoint(transfo, obj.Points);
            indices = round(points2(inds, :))+1;
            
            for i=1:length(inds)
                iInd = inds(i);
                
                % calcule jacobien pour points valides (repere image fixe)
                jac = parametricJacobian(transfo, obj.Points(iInd, :));
                
                % local gradient in moving image
                subs = num2cell(indices(i, :));
                grad = getPixel(obj.GradientImage, subs{:});
                
                % local contribution to metric gradient
                g(iInd,:) = grad*jac;
            end
            
            % compute gradient vectors weighted by local differences
            gd = g(inds,:).*diff(:, ones(1, nParams));
            
            % mean of valid gradient vectors
            grad = mean(gd, 1);
            
        end
            
        function [res, grad, isInside] = computeValueAndGradientLocal2d(obj)
            % Assumes gradient image is 2D.
            
            
            % error checking
            if isempty(obj.Transform)
                error('Gradient computation requires transform');
            end
            if isempty(obj.GradientImage)
                error('Gradient computation requires a gradient image');
            end
            
            % compute values in image 1
            [values1, inside1] = evaluate(obj.FixedImage, obj.Points);
            
            % compute values in image 2
            [values2, inside2] = evaluate(obj.TransformedImage, obj.Points);
            
            % keep only valid values
            isInside = inside1 & inside2;
            
            % compute result
            diff = values2(isInside) - values1(isInside);
            res = mean(diff.^2);
            
            %fprintf('Initial SSD: %f\n', res);
            
            % convert to indices
            inds = find(isInside);
            nbInds = length(inds);
            
            transfo = obj.Transform;
            nParams = length(getParameters(transfo));
            g = zeros(nbInds, nParams);
            
            % compute transformed coordinates
            points2 = transformPoint(transfo, obj.Points);
            
            % convert from physical coordinates to index coordinates
            % (assumes spacing is 1 and origin is 0)
            indices = round(points2(inds, :)) + 1;
            
            gradImg = obj.GradientImage.Data;
            
            for i = 1:length(inds)
                iInd = inds(i);
                
                % calcule jacobien pour points valides (repere image fixe)
                p0 = obj.Points(iInd, :);
                jac = parametricJacobian(transfo, p0);
                
                % local gradient in moving image
                ind1 = indices(i,1);
                ind2 = indices(i,2);
                
                grad = [gradImg(ind1, ind2, 1) gradImg(ind1, ind2, 2)];
                
                % local contribution to metric gradient
                g(iInd,:) = grad*jac;
            end
            
            % compute gradient vectors weighted by local differences
            gd = g(inds,:).*diff(:, ones(1, nParams));
            
            % mean of valid gradient vectors
            grad = mean(gd, 1);
        end
        
        
        function [res, grad, isInside] = computeValueAndGradientLocal3d(obj)
            % Assumes gradient image is 3D.
            
            
            % error checking
            if isempty(obj.Transform)
                error('Gradient computation requires transform');
            end
            if isempty(obj.GradientImage)
                error('Gradient computation requires a gradient image');
            end
            
            % compute values in image 1
            [values1, inside1] = evaluate(obj.FixedImage, obj.Points);
            
            % compute values in image 2
            [values2, inside2] = evaluate(obj.TransformedImage, obj.Points);
            
            % keep only valid values
            isInside = inside1 & inside2;
            
            % compute result
            diff = values2(isInside) - values1(isInside);
            res = mean(diff.^2);
            
            %fprintf('Initial SSD: %f\n', res);
            
            % convert to indices
            inds = find(isInside);
            nbInds = length(inds);
            
            transfo = obj.Transform;
            nParams = length(getParameters(transfo));
            g = zeros(nbInds, nParams);
            
            % compute transformed coordinates
            points2 = transformPoint(transfo, obj.Points);
            
            % convert from physical coordinates to index coordinates
            % (assumes spacing is 1 and origin is 0)
            indices = round(points2(inds, :)) + 1;
            
            gradImg = obj.GradientImage.Data;
            
            for i = 1:length(inds)
                iInd = inds(i);
                
                % calcule jacobien pour points valides (repere image fixe)
                p0 = obj.Points(iInd, :);
                jac = parametricJacobian(transfo, p0);
                
                % local gradient in moving image
                ind1 = indices(i,1);
                ind2 = indices(i,2);
                ind3 = indices(i,3);
                
                grad = [gradImg(ind1, ind2, ind3, 1) gradImg(ind1, ind2, ind3, 2) ...
                    gradImg(ind1, ind2, ind3, 3)];
                
                % local contribution to metric gradient
                g(iInd,:) = grad*jac;
            end
            
            % compute gradient vectors weighted by local differences
            gd = g(inds,:).*diff(:, ones(1, nParams));
            
            % mean of valid gradient vectors
            grad = mean(gd, 1);
            
        end
    end
end

end % classdef
