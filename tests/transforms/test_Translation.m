function test_suite = test_Translation(varargin)
%TEST_TRANSLATION  Test function for class Translation
%   output = test_Translation(input)
%
%   Example
%   test_Translation
%
%   See also
%
%
% ------
% Author: David Legland
% e-mail: david.legland@grignon.inra.fr
% Created: 2010-06-03,    using Matlab 7.9.0.529 (R2009b)
% Copyright 2010 INRA - Cepia Software Platform.

test_suite = buildFunctionHandleTestSuite(localfunctions);

function testEmptyConstructor %#ok<*DEFNU>

% test empty constructor
trans = Translation();
assertTrue(isvalid(trans));

function testConstructor2D

% test constructor with separate arguments
trans = Translation(2, 3);
assertTrue(isvalid(trans));

% test constructor with bundled arguments
trans = Translation([2, 3]);
assertTrue(isvalid(trans));

% test copy constructor
trans2 = Translation(trans);
assertTrue(isvalid(trans2));

function testConstructor3D

% test constructor with separate arguments
trans = Translation(2, 3, 4);
assertTrue(isvalid(trans));

% test constructor with bundled arguments
trans = Translation([2, 3, 4]);
assertTrue(isvalid(trans));

% test copy constructor
trans2 = Translation(trans);
assertTrue(isvalid(trans2));


function testIsa

trans = Translation([2 3]);
assertTrue(isa(trans, 'Transform'));
assertTrue(isa(trans, 'AffineTransform'));


function test_ToStruct
% Test call of function without argument

transfo = Translation([2 3]);
str = toStruct(transfo);
transfo2 = Translation.fromStruct(str);

assertTrue(isa(transfo2, 'Translation'));
assertElementsAlmostEqual(transfo2.U, transfo.U, 'absolute', .01);


function test_readWrite
% Test call of function without argument

% prepare
fileName = 'translation.transfo';
if exist(fileName, 'file')
    delete(fileName);
end

% arrange
transfo = Translation([2 3]);

% act
write(transfo, fileName);
transfo2 = Transform.read(fileName);

% assert
assertTrue(isa(transfo2, 'Translation'));
assertElementsAlmostEqual(transfo2.U, transfo.U, 'absolute', .01);

% clean up
delete(fileName);


