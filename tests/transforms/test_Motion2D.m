function test_suite = test_Motion2D
%TEST_MOTION2D  Test case for the file Motion2D
%
%   Test case for the file Motion2D
%
%   Example
%   test_Motion2D
%
%   See also
%   Motion2D
%
% ------
% Author: David Legland
% e-mail: david.legland@inra.fr
% Created: 2018-09-27,    using Matlab 9.4.0.813654 (R2018a)
% Copyright 2018 INRA - Cepia Software Platform.

test_suite = buildFunctionHandleTestSuite(localfunctions);

function test_Creation(testCase) %#ok<*DEFNU>
% Test call of function without argument
transfo = Motion2D(30, [10 20]);
assertTrue(isa(transfo, 'Motion2D'));

function test_ToStruct(testCase) %#ok<*DEFNU>
% Test call of function without argument

transfo = Motion2D(30, [10 20]);
str = toStruct(transfo);
transfo2 = Motion2D.fromStruct(str);

assertTrue(isa(transfo2, 'Motion2D'));
assertElementsAlmostEqual(transfo2.Theta, transfo.Theta, 'absolute', .01);
assertElementsAlmostEqual(transfo2.Translation, transfo.Translation, 'absolute', .01);


function test_readWrite(testCase) %#ok<*DEFNU>
% Test call of function without argument

% prepare
fileName = 'motion2D.transfo';
if exist(fileName, 'file')
    delete(fileName);
end

% arrange
transfo = Motion2D(30, [10 20]);

% act
write(transfo, fileName);
transfo2 = Transform.read(fileName);

% assert
assertTrue(isa(transfo2, 'Motion2D'));
assertElementsAlmostEqual(transfo2.Theta, transfo.Theta, 'absolute', .01);
assertElementsAlmostEqual(transfo2.Translation, transfo.Translation, 'absolute', .01);

% clean up
delete(fileName);

