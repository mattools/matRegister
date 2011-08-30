function test_suite = test_CenteredMotionTransform2D(varargin) %#ok<STOUT>
%TEST_CENTEREDMOTIONTRANSFORM2D  Test file for class CenteredMotionTransform2D
%   output = test_CenteredMotionTransform2D(input)
%
%   Example
%   test_CenteredMotionTransform2D
%
%   See also
%
%
% ------
% Author: David Legland
% e-mail: david.legland@grignon.inra.fr
% Created: 2010-06-17,    using Matlab 7.9.0.529 (R2009b)
% Copyright 2010 INRA - Cepia Software Platform.


initTestSuite;

function test_getAffineMatrix %#ok<*DEFNU>

center = [6 8];
T1 = Translation(-center);
rotMat = createRotation(deg2rad(30));
R = MatrixAffineTransform(rotMat);
T2 = Translation(center);

res = T2 * R * T1;

T = CenteredMotionTransform2D([30 0 0], 'center', center);

matTh = res.getAffineMatrix();
mat = T.getAffineMatrix();
assertElementsAlmostEqual(matTh, mat);


function test_getDimension

center = [6 8];
T = CenteredMotionTransform2D([30 0 0], 'center', center);

dim = getDimension(T);

assertEqual(2, dim);


function test_readWrite

% prepare
fileName = 'transfoFile.txt';
if exist(fileName, 'file')
    delete(fileName);
end

% create transfo
center = [6 8];
T = CenteredMotionTransform2D([30 0 0], 'center', center);
mat0 = getAffineMatrix(T);

% save the transfo
writeToFile(T, fileName);

% read a new transfo
T2 = CenteredMotionTransform2D.readFromFile(fileName);
mat2 = getAffineMatrix(T2);

assertElementsAlmostEqual(mat0, mat2);

% clean up
delete(fileName);

