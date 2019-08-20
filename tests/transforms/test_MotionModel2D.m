function testSuite = test_MotionModel2D(varargin)
%TEST_MOTIONMODEL2D  Test case for the file MotionModel2D
%
%   Test case for the file MotionModel2D

%   Example
%   test_MotionModel2D
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
MotionModel2D([0 0 0]);



function test_ToStruct
% Test call of function without argument

transfo = MotionModel2D([30 20 10]);
str = toStruct(transfo);
transfo2 = MotionModel2D.fromStruct(str);

assertTrue(isa(transfo2, 'MotionModel2D'));
assertElementsAlmostEqual(transfo2.Params, transfo.Params, 'absolute', .01);


function test_readWrite
% Test call of function without argument

% prepare
fileName = 'motionModel2D.transfo';
if exist(fileName, 'file')
    delete(fileName);
end

% arrange
transfo = MotionModel2D([30 20 10]);

% act
write(transfo, fileName);
transfo2 = Transform.read(fileName);

% assert
assertTrue(isa(transfo2, 'MotionModel2D'));
assertElementsAlmostEqual(transfo2.Params, transfo.Params, 'absolute', .01);

% clean up
delete(fileName);


