function test_suite = test_CenteredEulerTransform3D(varargin)
%test_CenteredEulerTransform3D  Test file for class CenteredMotionTransform2D
%   output = test_CenteredEulerTransform3D(input)
%
%   Example
%   test_CenteredEulerTransform3D
%
%   See also
%
%
% ------
% Author: David Legland
% e-mail: david.legland@inra.fr
% Created: 2010-06-17,    using Matlab 7.9.0.529 (R2009b)
% Copyright 2010 INRA - Cepia Software Platform.

test_suite = buildFunctionHandleTestSuite(localfunctions);

function test_getAffineMatrix %#ok<*DEFNU>

center = [6 8 9];
T1  = createTranslation3d(-center);
rotX = createRotationOx(deg2rad(3));
rotY = createRotationOy(deg2rad(4));
rotZ = createRotationOz(deg2rad(5));
R   = composeTransforms3d(rotX, rotY, rotZ);
T2  = createTranslation3d(center);
T0  = createTranslation3d([6 7 8]);
matTh = T0 * T2 * R * T1;

T = CenteredEulerTransform3D([3 4 5 6 7 8], 'center', center);

mat = T.getAffineMatrix();
assertElementsAlmostEqual(matTh, mat, 'absolute', .1);


function test_getDimension

center = [6 8 9];
T = CenteredEulerTransform3D([3 4 5 6 7 8], 'center', center);

dim = getDimension(T);

assertEqual(3, dim);


function test_readWrite

% prepare
fileName = 'transfoFile.txt';
if exist(fileName, 'file')
    delete(fileName);
end

% create transfo
center = [6 8 9];
T = CenteredEulerTransform3D([3 4 5 6 7 8], 'center', center);
mat0 = getAffineMatrix(T);

% save the transfo
writeToFile(T, fileName);

% read a new transfo
T2 = CenteredEulerTransform3D.readFromFile(fileName);
mat2 = getAffineMatrix(T2);

assertElementsAlmostEqual(mat0, mat2, 'absolute', .1);

% clean up
delete(fileName);

