function testSuite = test_SimilariryModel2D(varargin)
%TEST_SIMILARIRYMODEL2D  Test case for the file SimilariryModel2D
%
%   Test case for the file SimilariryModel2D

%   Example
%   test_SimilariryModel2D
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
SimilarityModel2D();



function test_ToStruct
% Test call of function without argument

transfo = SimilarityModel2D([30 20 10 -0.4]);
str = toStruct(transfo);
transfo2 = SimilarityModel2D.fromStruct(str);

assertTrue(isa(transfo2, 'SimilarityModel2D'));
assertElementsAlmostEqual(transfo2.params, transfo.params, 'absolute', .01);


function test_readWrite
% Test call of function without argument

% prepare
fileName = 'SimilarityModel2D.transfo';
if exist(fileName, 'file')
    delete(fileName);
end

% arrange
transfo = SimilarityModel2D([30 20 10 -0.4]);

% act
write(transfo, fileName);
transfo2 = Transform.read(fileName);

% assert
assertTrue(isa(transfo2, 'SimilarityModel2D'));
assertElementsAlmostEqual(transfo2.params, transfo.params, 'absolute', .01);

% clean up
delete(fileName);

