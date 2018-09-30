function testSuite = test_AffineTransformModel2D(varargin)
%TEST_AFFINETRANSFORMMODEL2D  Test case for the file AffineTransformModel2D
%
%   Test case for the file AffineTransformModel2D

%   Example
%   test_AffineTransformModel2D
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
AffineTransformModel2D();




function test_ToStruct
% Test call of function without argument

transfo = AffineTransformModel2D([1.1 0.2 0.3   0.1  1.2 0.5]);
str = toStruct(transfo);
transfo2 = AffineTransformModel2D.fromStruct(str);

assertTrue(isa(transfo2, 'AffineTransformModel2D'));
assertElementsAlmostEqual(transfo2.params, transfo.params, 'absolute', .01);


function test_readWrite
% Test call of function without argument

% prepare
fileName = 'AffineTransformModel2D.transfo';
if exist(fileName, 'file')
    delete(fileName);
end

% arrange
transfo = AffineTransformModel2D([1.1 0.2 0.3   0.1  1.2 0.5]);

% act
write(transfo, fileName);
transfo2 = Transform.read(fileName);

% assert
assertTrue(isa(transfo2, 'AffineTransformModel2D'));
assertElementsAlmostEqual(transfo2.params, transfo.params, 'absolute', .01);

% clean up
delete(fileName);

