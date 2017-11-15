function test_suite = testLinearInterpolator2D
% Test function for class LinearInterpolator2D
%   output = testImArea(input)
%
%   Example
%   testImEuler2d
%
%   See also
%

% ------
% Author: David Legland
% e-mail: david.legland@inra.fr
% Created: 2009-04-22,    using Matlab 7.7.0.471 (R2008b)
% Copyright 2009 INRA - Cepia Software Platform.

test_suite = buildFunctionHandleTestSuite(localfunctions);

function testCreateInterpolator %#ok<*DEFNU>

dim = [15 10];
REF = 200;

% create test image
dat = zeros([dim(2) dim(1)], 'uint8');
dat(3:end-3, 5:end-5) = REF;
img = Image.create(dat);

% create interpolator
interp = LinearInterpolator2D(img);
assertTrue(interp.isvalid);

% copy constructor
interp2 = LinearInterpolator2D(interp);
assertTrue(interp2.isvalid);

delete(interp);
assertTrue(interp2.isvalid);

function testEvaluateSingleCoord

REF = 200;
% evaluate for a single point
interp = createInterpolator();
x = 8; y = 5;
val = interp.evaluate(x, y);
assertElementsAlmostEqual(REF, val, 'absolute', .1);

function testEvaluateSinglePoint

REF = 200;
interp = createInterpolator();
x = 8; y = 5;
point = [x y];
val = interp.evaluate(point);
assertElementsAlmostEqual(REF, val, 'absolute', .1);

function testEvaluateCoordArray

lx = [6 7 8 9];
ly = [4 5];
[x, y] = meshgrid(lx, ly);

interp = createInterpolator();
val = interp.evaluate(x, y);
assertEqual(size(x), size(val));

[val, inside] = interp.evaluate(x, y);
assertElementsAlmostEqual(size(x), size(val));
assertElementsAlmostEqual(size(x), size(inside));

function testEvaluatePointArray

lx = [6 7 8 9];
ly = [4 5];
[x, y] = meshgrid(lx, ly);
points = [x(:) y(:)];
n = size(points, 1);

interp = createInterpolator();
[val, inside] = interp.evaluate(points);

assertElementsAlmostEqual(n, length(val));
assertElementsAlmostEqual(n, length(inside));


function testIsAnImageInterpolator

interp = createInterpolator();
assertTrue(isa(interp, 'ImageInterpolator2D'));

%% Utilitary functions

function interp = createInterpolator
dim = [15 10];
REF = 200;

% create test image
dat = zeros([dim(2) dim(1)], 'uint8');
dat(3:end-3, 5:end-5) = REF;
img = Image.create(dat);

% create interpolator
interp = LinearInterpolator2D(img);



