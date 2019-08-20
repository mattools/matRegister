%DEMORADIALSCALINGTRANSFORM2D  One-line description here, please.
%
%   output = demoRadialScalingTransform2D(input)
%
%   Example
%   demoRadialScalingTransform2D
%
%   See also
%
 
% ------
% Author: David Legland
% e-mail: david.legland@inra.fr
% Created: 2018-11-12,    using Matlab 8.6.0.267246 (R2015b)
% Copyright 2018 INRA - Cepia Software Platform.

% generate transform
angles = 0:359;
scalings = ones(size(angles));
scalings(angles > 45 & angles < 135) = 1.5;
scalings(angles > 225 & angles < 315) = 0.8;
transfo = RadialScalingTransform2D(angles, scalings);

% draw sample polygon
poly = resamplePolygon(rectToPolygon([-8 -5 16 10]), 1000);
figure; hold on; 
axis equal; axis([-10 10 -10 10]);
drawPolygon(poly, 'k');

% draw transformed polygon
poly2 = transformPoint(transfo, poly);
drawPolygon(poly2, 'm')