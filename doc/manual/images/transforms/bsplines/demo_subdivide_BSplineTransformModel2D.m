%DEMO_SUBDIVIDE_BSPLINETRANSFORMMODEL2D  One-line description here, please.
%
%   output = demo_subdivide_BSplineTransformModel2D(input)
%
%   Example
%   demo_subdivide_BSplineTransformModel2D
%
%   See also
%
 
% ------
% Author: David Legland
% e-mail: david.legland@inrae.fr
% INRAE - BIA Research Unit - BIBS Platform (Nantes)
% Created: 2026-03-23,    using Matlab 25.1.0.2973910 (R2025a) Update 1
% Copyright 2026 INRAE.

% create a new transform, based on 3-by-3 grid, with spacing 60 and origin
% at point (40,40).
transfo = BSplineTransformModel2D([3 3], [60 60], [40 40]);

transfo.Params = [...
    +1 +1   -1 +1    0  0 ...
    +1 -1   -1 -1   +1 -1 ...
     0  0   -1 +1   +1 +1 ...
    ] * 25;

figure; hold on; axis equal; axis([0 200 0 200]);
drawGrid(transfo);
drawVertexShifts(transfo, 'b');

print(gcf, 'BSplineTransformModel2D_3x3.png', '-dpng');

% compute a subdivide transform, 
% with new grid spacingf equal to 40
transfo2 = subdivide(transfo);

figure; hold on; axis equal; axis([0 200 0 200]);
drawGrid(transfo2);
drawVertexShifts(transfo2, 'b');

print(gcf, 'BSplineTransformModel2D_3x3_subdivided.png', '-dpng');
