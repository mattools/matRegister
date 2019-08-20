function testSuite = test_RadialScalingTransform2D(varargin)
%TEST_RADIALSCALINGTRANSFORM2D  Test case for the file RadialScalingTransform2D
%
%   Test case for the file RadialScalingTransform2D

%   Example
%   test_RadialScalingTransform2D
%
%   See also
%

% ------
% Author: David Legland
% e-mail: david.legland@inra.fr
% Created: 2018-11-12,    using Matlab 8.6.0.267246 (R2015b)
% Copyright 2018 INRA - Cepia Software Platform.

testSuite = buildFunctionHandleTestSuite(localfunctions);

function test_EmptyConstructor %#ok<*DEFNU>
% Test call of function without argument
RadialScalingTransform2D();


function test_ToStruct
% Test call of function without argument

angles = 0:5:359;
scalings = ones(size(angles)) + 0.2 * cos(deg2rad(2*angles));
transfo = RadialScalingTransform2D(angles, scalings);
str = toStruct(transfo);
transfo2 = RadialScalingTransform2D.fromStruct(str);

assertTrue(isa(transfo2, 'RadialScalingTransform2D'));
assertElementsAlmostEqual(transfo2.Angles, transfo.Angles, 'absolute', .01);
assertElementsAlmostEqual(transfo2.Scalings, transfo.Scalings, 'absolute', .01);


function test_readWrite
% Test call of function without argument

% prepare
fileName = 'RadialScalingTransform2D.transfo';
if exist(fileName, 'file')
    delete(fileName);
end

% arrange
angles = 0:5:359;
scalings = ones(size(angles)) + 0.2 * cos(deg2rad(2*angles));
transfo = RadialScalingTransform2D(angles, scalings);

% act
write(transfo, fileName);
transfo2 = Transform.read(fileName);

% assert
assertTrue(isa(transfo2, 'RadialScalingTransform2D'));
assertElementsAlmostEqual(transfo2.Angles, transfo.Angles, 'absolute', .01);
assertElementsAlmostEqual(transfo2.Scalings, transfo.Scalings, 'absolute', .01);

% clean up
delete(fileName);



