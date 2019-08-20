function [res, grad, isInside] = computeValueAndGradient(obj, varargin)
% Compute metric value and gradient
%
%   [RES, DERIV] = obj.computeValueAndGradient();
%   This syntax requires that both fields 'transform' and 'gradientImage'
%   have been initialized.
%
%   [RES, DERIV] = obj.computeValueAndGradient(TRANSFO, GRADX, GRADY);
%   [RES, DERIV] = obj.computeValueAndGradient(TRANSFO, GRADX, GRADY, GRADZ);
%   This (deprecated) syntax passes transform model and gradient components
%   as input arguments.
%
% Example:
% transfo = Translation2DModel([1.2 2.3]);
% ssdMetric = SumOfSquaredDifferencesMetric(img1, img2, transfo);
% res = ssdMetric.computeValueAndGradient(model);
%

nd = ndims(obj.Img1);

% The first part of the file consists in analyzing input, and to call the
% most appropriate subfunction
if ~isempty(obj.Transform) && ~isempty(obj.GradientImage)
    
    % If the gradient image is an image function, it is assumed to be a 
    % gradient interpolator or a gradient evaluator
    if isa(obj.GradientImage, 'ImageFunction')
        [res, grad, isInside] = ...
            computeValueAndGradientFromGradientFunction(obj);
        return;
    end
    
    % if gradient image is a standard image, use specific methods that
    % perform dimension-specific nearest-neighbor interpolation
    if nd == 2
        [res, grad, isInside] = computeValueAndGradientLocal2d(obj);
    elseif nd == 3
        [res, grad, isInside] = computeValueAndGradientLocal3d(obj);
    else
        [res, grad, isInside] = computeValueAndGradientLocal(obj);
    end
    
else
    % deprecation warning
    warning('matRegister:deprecated', ...
        'Deprecated syntax. Please initialize metric fields instead');
  
    if length(varargin) < nd
        error('Requires as many gradient components as the number of dimensions');
    end
    
    % assumes transform and gradient components are given as arguments
    if nd == 2
        [res, grad, isInside] = computeValueAndGradient2d(obj, varargin{:});
    else
        [res, grad, isInside] = computeValueAndGradient3d(obj, varargin{:});
    end
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
[values1, inside1] = evaluate(obj.Img1, obj.Points);

% compute values in image 2
[values2, inside2] = evaluate(obj.Img2, obj.Points);

% keep only valid values
isInside = inside1 & inside2;

% compute result
diff = values2(isInside) - values1(isInside);

% average over all points
np = length(isInside);
res = sum(diff .^ 2) / np;

%fprintf('Initial SSD: %f\n', res);

% convert to indices
inds    = find(isInside);
nInds   = length(inds);

transfo = obj.Transform;
nParams = length(getParameterLength(transfo));
gd = zeros(nInds, nParams);

% convert from physical coordinates to index coordinates
% (assumes spacing is 1 and origin is 0)
points2 = transformPoint(transfo, obj.Points);
indices = round(points2(inds, :)) + 1;

for i = 1:length(inds)
    iInd = inds(i);
    
    % compute jacobian for valid points (in fixed image reference system)
    jac = parametricJacobian(transfo, obj.Points(iInd, :));
    
    % local gradient in moving image
    subs = num2cell(indices(i, :));
    grad = getPixel(obj.GradientImage, subs{:});
    
    % local contribution to metric gradient
    gd(i, :) = grad * jac;
end

% compute gradient vectors weighted by local differences
gd = gd .* diff(:, ones(1, nParams));

% mean of valid gradient vectors
grad = mean(gd, 1);





function [res, grad, isInside] = computeValueAndGradientLocal2d(obj)
% Compute metric value and gradient in 2D, using inner gradient image
%
% Assumes that gradient is a 2D image.

% error checking
if isempty(obj.Transform)
    error('Gradient computation requires transform');
end
if isempty(obj.GradientImage)
    error('Gradient computation requires a gradient image');
end

% compute values in image 1
[values1, inside1] = evaluate(obj.Img1, obj.Points);

% compute values in image 2
[values2, inside2] = evaluate(obj.Img2, obj.Points);

% keep only valid values
isInside = inside1 & inside2;

% compute result
diff = values2 - values1;

% average over all points
% np  = length(isInside);
% res = sum(diff .^ 2) / np;
res = mean(diff .^ 2);

%fprintf('Initial MSD: %f\n', res);

% convert to indices
inds    = find(isInside);
nInds   = length(inds);

transfo = obj.Transform;
nParams = getParameterLength(transfo);
gd      = zeros(nInds, nParams);

% compute transformed coordinates
points2 = transformPoint(transfo, obj.Points);

% convert from physical coordinates to index coordinates
% (assumes spacing is 1 and origin is 0)
% indices = round(points2(inds, :)) + 1;
indices = pointToIndex(obj.GradientImage, points2(inds,:));

gradImg = obj.GradientImage.Data;

% iterate over sampling points within bounds
for i = 1:nInds
    iInd = inds(i);
    
    % compute jacobian for valid points (in fixed image reference system)
    p0 = obj.Points(iInd, :);
    jac = parametricJacobian(transfo, p0);
    
    % local gradient in moving image
    ind1 = indices(i,1);
    ind2 = indices(i,2);
    
    grad = [gradImg(ind1, ind2, 1, 1) gradImg(ind1, ind2, 1, 2)];
    
    % local contribution to metric gradient
    gd(i, :) = grad * jac;
end

% compute gradient vectors weighted by local differences
% gd = gd .* diff(inds, ones(1, nParams));
gd = bsxfun(@times, gd, diff(inds));

% mean of valid gradient vectors
grad = mean(gd, 1);



function [res, grad, isInside] = computeValueAndGradientLocal3d(obj)
%Assumes gradient image is 3D

% TODO: need to update to compute diff on all values

% compute values in image 1
[values1, inside1] = evaluate(obj.Img1, obj.Points);

% compute values in image 2
[values2, inside2] = evaluate(obj.Img2, obj.Points);

% keep only valid values
isInside = inside1 & inside2;

if sum(isInside) < 100
    error('Too many points outside registration window');
end

% compute result
diff = values2(isInside) - values1(isInside);

% average over all points
np  = length(isInside);
res = sum(diff.^2)/np;

%fprintf('Initial SSD: %f\n', res);

% convert to indices
inds    = find(isInside);
nInds   = length(inds);

transfo = obj.Transform;
nParams = getParameterLength(transfo);
gd = zeros(nInds, nParams);

% compute transformed coordinates
points2 = transformPoint(transfo, obj.Points);

% convert from physical coordinates to index coordinates
% (assumes spacing is 1 and origin is 0)
% indices = round(points2(inds, :)) + 1;
gradImg = obj.GradientImage.Data;
indices = pointToIndex(gradImg, points2(inds,:));

% iterate over sampling points within bounds
for i = 1:length(inds)
    iInd = inds(i);
    
    % calcule jacobien pour points valides (repere image fixe)
    p0 = obj.Points(iInd, :);
    jac = parametricJacobian(transfo, p0);
    
    % local gradient in moving image
    ind1 = indices(i,1);
    ind2 = indices(i,2);
    ind3 = indices(i,3);

    grad = [...
        gradImg(ind1, ind2, ind3, 1) ...
        gradImg(ind1, ind2, ind3, 2) ...
        gradImg(ind1, ind2, ind3, 3) ];

    % local contribution to metric gradient
    tmp = grad * jac;
    gd(i, :) = tmp;
end

% compute gradient vectors weighted by local differences
gd = gd .* diff(:, ones(1, nParams));

% mean of valid gradient vectors
grad = mean(gd, 1);


function [res, grad, isInside] = computeValueAndGradientFromGradientFunction(obj)
% Assumes gradient image is an ImageFunction that evaluates or interpolates
% gradient of another image
% the interpolated/evaluated gradient is not transformed, obj operation is
% left to obj function

% compute values in image 1
[values1, inside1] = evaluate(obj.Img1, obj.Points);

% compute values in image 2
[values2, inside2] = evaluate(obj.Img2, obj.Points);

% keep only valid values
isInside = inside1 & inside2;

if sum(isInside) < 100
    error('Too many points outside registration window');
end

% compute result
diff = values2(isInside) - values1(isInside);

% average over all points
np  = length(isInside);
res = sum(diff .^ 2) / np;

%fprintf('Initial SSD: %f\n', res);


transfo = obj.Transform;
nParams = getParameterLength(transfo);

% compute transformed coordinates
points2 = transformPoint(transfo, obj.Points);

% evaluate gradient, and re-compute points within image frame, as gradient
% evaluator can have different behaviour at image borders.
[gradVals, gradInside] = evaluate(obj.GradientImage, points2);

% convert to indices
inds    = find(gradInside);
nInds   = length(inds);
gd      = zeros(nInds, nParams);

for i = 1:nInds
    iInd = inds(i);
    
    % compute jacobian for valid points (in fixed image reference system)
    p0  = obj.Points(iInd, :);
    jac = parametricJacobian(transfo, p0);
    
    % % local contribution to metric gradient
    gd(i, :) = gradVals(iInd, :)*jac;
end

% re-compute differences, by considering position that can be used for
% computing gradient
diff = values2(gradInside) - values1(gradInside);

% compute gradient vectors weighted by local differences
gd = gd .* diff(:, ones(1, nParams));

% remove some NAN values that could occur for an obscure reason
gd = gd(~isnan(gd(:,1)), :);

% mean of valid gradient vectors
grad = mean(gd, 1);



function [res, grad, isInside] = computeValueAndGradient2d(obj, transfo, gx, gy)
% Old function to compute metric and gradient of a 2D image using 2 args

% compute values in image 1
[values1, inside1] = evaluate(obj.Img1, obj.Points);

% compute values in image 2
[values2, inside2] = evaluate(obj.Img2, obj.Points);

% keep only valid values
isInside = inside1 & inside2;

% compute result
diff = values2(isInside) - values1(isInside);

% average over all points
np = length(isInside);
res = sum(diff .^ 2) / np;

%fprintf('Initial SSD: %f\n', res);


% convert to indices
inds    = find(isInside);
nInds   = length(inds);

%nPoints = size(points, 1);
nParams = getParameterLength(transfo);
gd = zeros(nInds, nParams);

% convert from physical coordinates to index coordinates
% (assumes spacing is 1 and origin is 0)
% also converts from (x,y) to (i,j)
points2 = TransformPoint(transfo, obj.Points);
index = round(points2(inds, [2 1]))+1;

for i = 1:nInds
    % compute jacobian for valid points (in fixed image reference system)
    jac = parametricJacobian(transfo, obj.Points(inds(i), :));
    
    % local gradient in moving image
    i1 = index(i, 1);
    i2 = index(i, 2);
    grad = [gx(i1, i2) gy(i1, i2)];
    
    % local contribution to metric gradient
    gd(i, :) = grad * jac;
end

% calcul du vecteur gradient pondere par difference locale
gd = gd .* diff(:, ones(1, nParams));

% somme des vecteurs gradient valides
grad = mean(gd, 1);


function [res, grad, isInside] = computeValueAndGradient3d(obj, transfo, gx, gy, gz)
% Old function to compute metric and gradient of a 3D image using 3 args

% compute values in image 1
[values1, inside1] = evaluate(obj.Img1, obj.Points);

% compute values in image 2
[values2, inside2] = evaluate(obj.Img2, obj.Points);

% keep only valid values
isInside = inside1 & inside2;

% compute result
diff = values2(isInside) - values1(isInside);

% average over all points
np  = length(isInside);
res = sum(diff .^ 2) / np;


% convert to indices
inds    = find(isInside);
nInds   = length(inds);

nParams = getParameterLength(transfo);
gd      = zeros(nInds, nParams);

% convert from physical coordinates to index coordinates
% (assumes spacing is 1 and origin is 0)
% also converts from (x,y) to (i,j)
points2 = transformPoint(transfo, obj.Points);
index = round(points2(inds, [2 1 3]))+1;

for i = 1:nInds
    % compute jacobian for valid points (in fixed image reference system)
    jac = parametricJacobian(transfo, obj.Points(inds(i),:));
    
    % local gradient in moving image
    i1 = index(i, 1);
    i2 = index(i, 2);
    i3 = index(i, 3);
    grad = [gx(i1,i2,i3) gy(i1,i2,i3) gz(i1,i2,i3)];
    
    % local contribution to metric gradient
    gd(i, :) = grad * jac;
end

% calcul du vecteur gradient pondere par difference locale
gd = gd .* diff(:, ones(1, nParams));

% somme des vecteurs gradient valides
grad = mean(gd, 1);
