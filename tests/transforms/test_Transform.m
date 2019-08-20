function test_suite = test_Transform(varargin)
%testTransform  One-line description here, please.
%   output = testTransform(input)
%
%   Example
%   testTransform
%
%   See also
%
%
% ------
% Author: David Legland
% e-mail: david.legland@grignon.inra.fr
% Created: 2010-06-17,    using Matlab 7.9.0.529 (R2009b)
% Copyright 2010 INRA - Cepia Software Platform.

test_suite = buildFunctionHandleTestSuite(localfunctions);

function test_compose %#ok<*DEFNU>

% transform parameters
center = [6 8];
angle = deg2rad(30);

% create transform objects
T1 = Translation(-center);
rotMat = createRotation(angle);
R = MatrixAffineTransform(rotMat);
T2 = Translation(center);

% get transform matrices
t1Mat = affineMatrix(T1);
t2Mat = affineMatrix(T2);
resMat = t2Mat*rotMat*t1Mat;

% create composed transforms
res1 = MatrixAffineTransform(resMat);
res2 = T2.compose(R).compose(T1);
res3 = CenteredMotionTransform2D([30 0 0], 'Center', center);

%% transform a set of points using both transforms
pts0 = [5 6; 3 4;-1 2];
pts1 = transformPoint(res1, pts0);
pts2 = transformPoint(res2, pts0);
pts3 = transformPoint(res3, pts0);

% transformed points should be the same
assertElementsAlmostEqual(pts1, pts2, 'absolute', .1);
assertElementsAlmostEqual(pts1, pts3, 'absolute', .1);
