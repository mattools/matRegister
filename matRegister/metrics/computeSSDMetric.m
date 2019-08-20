function res = computeSSDMetric(img1, img2, points)
% Sum of squared differences between two (interpolated) images.
%   output = computeSSDMetric(input)
%
%   Example
%   computeSSDMetric
%
%   See also
%

% ------
% Author: David Legland
% e-mail: david.legland@inra.fr
% Created: 2010-06-10,    using Matlab 7.9.0.529 (R2009b)
% Copyright 2010 INRA - Cepia Software Platform.

if isa(img1, 'ImageFunction')
    % ok
elseif isa(img1, 'Image')
    % if points are no specified, compute their position from image
    if nargin==2
        x = getX(img1);
        y = getY(img1);
        points = [x(:) y(:)];
    end
    % convert image to interpolator
    img1 = LinearInterpolator2D(img1);
else
    error('First argument must be an image function');
end

if isa(img2, 'ImageFunction')
    % ok
elseif isa(img2, 'Image')
    % convert image to interpolator
    img2 = LinearInterpolator2D(img2);
else
    error('Second argument must be an image function');
end
 
% compute values in image 1
[values1, inside1] = evaluate(img1, points);

% compute values in image 2
[values2, inside2] = evaluate(img2, points);

% keep only valid values
inds = inside1 & inside2;

% compute result
diff = (values2(inds) - values1(inds)).^2;
res = sum(diff);
