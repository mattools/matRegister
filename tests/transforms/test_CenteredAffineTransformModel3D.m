function test_suite = test_CenteredAffineTransformModel3D(varargin)
%test_CenteredAffineTransformModel3D  One-line description here, please.
%   output = test_CenteredAffineTransformModel3D(input)
%
%   Example
%   test_CenteredAffineTransformModel3D
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

function test_createFromMatrix %#ok<*DEFNU>

mat = [diag([2 3 4]) [1;2;3] ; 0 0 0 1];
trans = CenteredAffineTransformModel3D(mat);
mat2 = trans.getAffineMatrix;

assertEqual(mat, mat2);

trans = CenteredAffineTransformModel3D(mat(1:3, :));
mat2 = trans.getAffineMatrix;

assertEqual(mat, mat2);


function test_createFromScalar 

trans = CenteredAffineTransformModel3D(3);
mat = trans.getAffineMatrix();
assertEqual(eye(4), mat);


function test_createCenteredTranslation

% create transfo
vect    = [3 4 5];
trans   = createTranslation3d(vect);
transfo = CenteredAffineTransformModel3D(trans);

% add center
center  = [10 20 30];
transfo.setCenter(center);

% center should not change translation
mat = transfo.getAffineMatrix();
assertEqual(trans, mat);


function test_createFromMotion

% create transfo
rot   = eulerAnglesToRotation3d(10, 20, 30);
crot  = recenterTransform3d(rot, [3 4 5]);
transfo = CenteredAffineTransformModel3D(crot);

% center should not change translation
mat = transfo.getAffineMatrix();
assertEqual(crot, mat);


function test_createCenteredRotation

% create transfo
rot   = eulerAnglesToRotation3d(10, 20, 30);
transfo = CenteredAffineTransformModel3D(rot);

% add center
center  = [10 20 30];
transfo.setCenter(center);

% center should not change translation
assertElementsAlmostEqual(center, transfo.transformPoint(center), 'absolute', .1);

