function test_suite = test_AffineTransform(varargin)
%test_AffineTransform  One-line description here, please.
%   output = test_AffineTransform(input)
%
%   Example
%   testAffineTransform
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

function test_createTranslation_2d %#ok<*DEFNU>

mat = affineMatrix(AffineTransform.createTranslation([2 3]));
matTh = [1 0 2 ;0 1 3; 0 0 1];
assertElementsAlmostEqual(matTh, mat, 'absolute', .1);

function test_createScaling_2d %#ok<*DEFNU>

mat = affineMatrix(AffineTransform.createScaling([3 2]));
matTh = [3 0 0 ;0 2 0; 0 0 1];
assertElementsAlmostEqual(matTh, mat, 'absolute', .1);

function test_mtimes

% Compose two translations
T1 = Translation([2 3]);
T2 = Translation([4 5]);

res = T1*T2;
mat = affineMatrix(res);

matTh = [1 0 6;0 1 8;0 0 1];
assertElementsAlmostEqual(matTh, mat, 'absolute', .1);


center = [6 8];
T1 = Translation(-center);
rotMat = createRotation(deg2rad(30));
R = MatrixAffineTransform(rotMat);
T2 = Translation(center);

res = T2*R*T1;

T = CenteredMotionTransform2D([30 0 0], 'center', center);

assertElementsAlmostEqual(affineMatrix(res), affineMatrix(T), 'absolute', .1);
