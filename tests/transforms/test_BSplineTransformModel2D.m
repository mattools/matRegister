function testSuite = test_BSplineTransformModel2D(varargin)
%TEST_BSPLINETRANSFORMMODEL2D  Test case for the file BSplineTransformModel2D
%
%   Test case for the file BSplineTransformModel2D

%   Example
%   test_BSplineTransformModel2D
%
%   See also
%

% ------
% Author: David Legland
% e-mail: david.legland@inra.fr
% Created: 2018-09-30,    using Matlab 8.6.0.267246 (R2015b)
% Copyright 2018 INRA - Cepia Software Platform.

testSuite = buildFunctionHandleTestSuite(localfunctions);

function test_Simple %#ok<*DEFNU>
% Test call of function without argument
BSplineTransformModel2D([3 3], [10 10], [0 0]);


function test_ToStruct
% Test call of function without argument

transfo = BSplineTransformModel2D([3 3], [10 10], [0 0]);
transfo.params = zeros(1, 18);
str = toStruct(transfo);
transfo2 = BSplineTransformModel2D.fromStruct(str);

assertTrue(isa(transfo2, 'BSplineTransformModel2D'));
assertElementsAlmostEqual(transfo2.params, transfo.params, 'absolute', .01);


function test_readWrite
% Test call of function without argument

% prepare
fileName = 'BSplineTransformModel2D.transfo';
if exist(fileName, 'file')
    delete(fileName);
end

% arrange
transfo = BSplineTransformModel2D([3 3], [10 10], [0 0]);
transfo.params = zeros(1, 18);

% act
write(transfo, fileName);
transfo2 = Transform.read(fileName);

% assert
assertTrue(isa(transfo2, 'BSplineTransformModel2D'));
assertElementsAlmostEqual(transfo2.params, transfo.params, 'absolute', .01);

% clean up
delete(fileName);



