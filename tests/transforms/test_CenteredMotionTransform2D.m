function test_suite = test_CenteredMotionTransform2D(varargin)
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

test_suite = buildFunctionHandleTestSuite(localfunctions);

function test_affineMatrix %#ok<*DEFNU>

center = [6 8];
T1 = Translation(-center);
rotMat = createRotation(deg2rad(30));
R = MatrixAffineTransform(rotMat);
T2 = Translation(center);

res = T2 * R * T1;

T = CenteredMotionTransform2D([30 0 0], 'center', center);

matTh = affineMatrix(res);
mat = affineMatrix(T);
assertElementsAlmostEqual(matTh, mat, 'absolute', .1);


function test_getDimension

center = [6 8];
T = CenteredMotionTransform2D([30 0 0], 'center', center);

dim = getDimension(T);

assertEqual(2, dim);


function test_writeToFile

% prepare
fileName = 'transfoFile.txt';
if exist(fileName, 'file')
    delete(fileName);
end

% create transfo
center = [6 8];
T = CenteredMotionTransform2D([30 0 0], 'center', center);
mat0 = affineMatrix(T);

% save the transfo
writeToFile(T, fileName);

% read a new transfo
T2 = CenteredMotionTransform2D.readFromFile(fileName);
mat2 = affineMatrix(T2);

assertElementsAlmostEqual(mat0, mat2, 'absolute', .1);

% clean up
delete(fileName);


function test_ToStruct
% Test call of function without argument

transfo = CenteredMotionTransform2D([10 20 30], 'center', [50 50]);
str = toStruct(transfo);
transfo2 = CenteredMotionTransform2D.fromStruct(str);

assertTrue(isa(transfo2, 'CenteredMotionTransform2D'));
assertElementsAlmostEqual(transfo2.params, transfo.params, 'absolute', .01);


function test_readWrite
% Test call of function without argument

% prepare
fileName = 'CenteredMotionTransform2D.transfo';
if exist(fileName, 'file')
    delete(fileName);
end

% arrange
transfo = CenteredMotionTransform2D([10 20 30], 'center', [50 50]);

% act
write(transfo, fileName);
transfo2 = Transform.read(fileName);

% assert
assertTrue(isa(transfo2, 'CenteredMotionTransform2D'));
assertElementsAlmostEqual(transfo2.params, transfo.params, 'absolute', .01);
assertElementsAlmostEqual(transfo2.center, transfo.center, 'absolute', .01);

% clean up
delete(fileName);

