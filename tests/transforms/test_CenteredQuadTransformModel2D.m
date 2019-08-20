function testSuite = test_CenteredQuadTransformModel2D(varargin)
%TEST_CENTEREDQUADTRANSFORMMODEL2D  Test case for the file CenteredQuadTransformModel2D
%
%   Test case for the file CenteredQuadTransformModel2D

%   Example
%   test_CenteredQuadTransformModel2D
%
%   See also
%

% ------
% Author: David Legland
% e-mail: david.legland@inra.fr
% Created: 2018-10-01,    using Matlab 8.6.0.267246 (R2015b)
% Copyright 2018 INRA - Cepia Software Platform.

testSuite = buildFunctionHandleTestSuite(localfunctions);

function test_Simple %#ok<*DEFNU>
% Test call of function without argument
CenteredQuadTransformModel2D();


function test_ToStruct
% Test call of function without argument

params = [.1 .2  1.1 0.2  0.1 1.2  0.01 0.01 0.01 0.01 0.01 0.01];
transfo = CenteredQuadTransformModel2D(params, 'Center', [50 50]);
str = toStruct(transfo);
transfo2 = CenteredQuadTransformModel2D.fromStruct(str);

assertTrue(isa(transfo2, 'CenteredQuadTransformModel2D'));
assertElementsAlmostEqual(transfo2.Params, transfo.Params, 'absolute', .01);


function test_readWrite
% Test call of function without argument

% prepare
fileName = 'CenteredQuadTransformModel2D.transfo';
if exist(fileName, 'file')
    delete(fileName);
end

% arrange
params = [.1 .2  1.1 0.2  0.1 1.2  0.01 0.01 0.01 0.01 0.01 0.01];
transfo = CenteredQuadTransformModel2D(params, 'Center', [50 50]);

% act
write(transfo, fileName);
transfo2 = Transform.read(fileName);

% assert
assertTrue(isa(transfo2, 'CenteredQuadTransformModel2D'));
assertElementsAlmostEqual(transfo2.Params, transfo.Params, 'absolute', .01);
assertElementsAlmostEqual(transfo2.Center, transfo.Center, 'absolute', .01);

% clean up
delete(fileName);

