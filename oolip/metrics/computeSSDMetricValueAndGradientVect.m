function [res grad] = computeSSDMetricValueAndGradientVect(img1, img2, points, transfo, grad)
%COMPUTESSDMETRIC Compute SSD value and gradient in 2D image using gradient image
%   output = computeSSDMetric(input)
%
%   Version utilisant une image gradient
%   Example
%   computeSSDMetric
%
%   See also
%
%
% ------
% Author: David Legland
% e-mail: david.legland@grignon.inra.fr
% Created: 2010-06-10,    using Matlab 7.9.0.529 (R2009b)
% Copyright 2010 INRA - Cepia Software Platform.


%% Process input arguments

if isa(img1, 'ImageFunction')
    % ok
elseif isa(img1, 'Image')
    % if points are no specified, compute their position from image
    if nargin == 2
        x = img1.getX();
        y = img1.getY();
        points = [x(:) y(:)];
    end
    % convert image to interpolator
    disp('convert img1 to interpolator');
    img1 = LinearInterpolator2D(img1);
else
    error('First argument must be an image function');
end

if isa(img2, 'ImageFunction')
    % ok
elseif isa(img2, 'Image')
    % convert image to interpolator
    disp('convert img2 to interpolator');
    img2 = LinearInterpolator2D(img2);
else
    error('Second argument must be an image function');
end


%% Compute metric value 

% compute values in image 1
[values1 inside1] = img1.evaluate(points);

% compute values in image 2
[values2 inside2] = img2.evaluate(points);

% keep only valid values
insideBoth = inside1 & inside2;

% compute result
diff = values2(insideBoth) - values1(insideBoth);
res = mean(diff .^ 2);


%% Compute gradient direction

% convert to indices
inds    = find(insideBoth);
nbInds  = length(inds);

nParams = length(transfo.getParameters());
g = zeros(nbInds, nParams);

% convert from physical coordinates to index coordinates
% (assumes spacing is 1 and origin is 0)
points2 = transfo.transformPoint(points);

% just need to keep integer coordinates
indices = round(points2(inds, :));

% iterate on valid points
for i = 1:length(inds)
    % calcule jacobien pour points valides (repere image fixe)
    jac = transfo.getParametricJacobian(points(inds(i), :));

    grad_i = grad.getPixel(indices(i,1), indices(i,2));
    
    % local contribution to metric gradient
    g(inds(i),:) = grad_i * jac;
end

% calcul du vecteur gradient pondere par difference locale
gd = g(inds,:) .* diff(:, ones(1, nParams));

% average of valid gradient vectors
grad = mean(gd, 1);

