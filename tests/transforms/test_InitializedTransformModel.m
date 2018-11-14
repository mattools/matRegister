function test_suite = test_InitializedTransformModel
%TEST_INITIALIZEDTRANSFORMMODEL  Test case for the file InitializedTransformModel
%
%   Test case for the file InitializedTransformModel
%
%   Example
%   test_InitializedTransformModel
%
%   See also
%   InitializedTransformModel

% ------
% Author: David Legland
% e-mail: david.legland@inra.fr
% Created: 2018-10-17,    using Matlab 9.5.0.944444 (R2018b)
% Copyright 2018 INRA - Cepia Software Platform.

test_suite = buildFunctionHandleTestSuite(localfunctions);

function test_Simple(testCase) %#ok<*DEFNU>
% Test call of function without argument

init = MatrixAffineTransform(eye(3, 3));
trans = BSplineTransformModel2D([3 3], [10 10], [0 0]);
transfo = InitializedTransformModel(init, trans); %#ok<NASGU>

function test_ReadWrite(testCase)
% Test call of function without argument

init = MatrixAffineTransform(eye(3, 3));
trans = BSplineTransformModel2D([3 3], [10 10], [0 0]);
transfo = InitializedTransformModel(init, trans);

transfo2 = InitializedTransformModel.fromStruct(toStruct(transfo));
assertTrue(isa(transfo2, 'InitializedTransformModel'));
assertTrue(isa(transfo2.initial, 'Transform'));
assertTrue(isa(transfo2.transform, 'ParametricTransform'));
