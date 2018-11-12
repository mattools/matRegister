%DEMO_RADIALSCALINGTRANSFORM2D_CAMERAMAN  One-line description here, please.
%
%   output = demo_RadialScalingTransform2D_cameraman(input)
%
%   Example
%   demo_RadialScalingTransform2D_cameraman
%
%   See also
%
 
% ------
% Author: David Legland
% e-mail: david.legland@inra.fr
% Created: 2018-11-12,    using Matlab 8.6.0.267246 (R2015b)
% Copyright 2018 INRA - Cepia Software Platform.

img = imread('cameraman.tif');

% generate transform
angles = 0:359;
scalings = 1 + 0.2 * cos(deg2rad(2 * angles));
transfo = RadialScalingTransform2D(angles, scalings);

figure(1); close;
figure(1); imshow(img);
hold on; 


% %% few points in the middle
% 
% lx = 100:5:155;
% ly = 100:5:155;
% [x, y] = meshgrid(lx, ly);
% pts = [x(:) y(:)];
% drawPoint(pts, 'b.');
% 
% pts2 = transformPoint(transfo, pts);
% drawPoint(pts2, 'go');


%% points within the whole image

center = [128 128];
lx = 10:5:250;
ly = 10:5:250;
[x, y] = meshgrid(lx, ly);
pts = [x(:)-center(1) y(:)-center(1)];

pts2 = transformPoint(transfo, pts);
pts2 = bsxfun(@plus, pts2, center);
drawPoint(pts2, 'g.');


%% reamsple image

% resample the whole image
lx = 1:size(img, 2);
ly = 1:size(img, 1);
[x, y] = meshgrid(lx, ly);
pts = [x(:)-center(1) y(:)-center(1)];

% transform points
pts2 = transformPoint(transfo, pts);
pts2 = bsxfun(@plus, pts2, center);

% create new image
img2 = reshape(imEvaluate(img, pts2), size(img));

% display result
figure;
imshow(img2, [0 255]);
