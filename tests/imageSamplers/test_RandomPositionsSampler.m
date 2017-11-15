function test_suite = test_RandomPositionsSampler(varargin)
%TEST_FULLIMAGESAMPLER  One-line description here, please.
%
%   output = test_RandomPositionsSampler(input)
%
%   Example
%   test_RandomPositionsSampler
%
%   See also
%
%
% ------
% Author: David Legland
% e-mail: david.legland@grignon.inra.fr
% Created: 2011-07-20,    using Matlab 7.9.0.529 (R2009b)
% Copyright 2011 INRA - Cepia Software Platform.

test_suite = buildFunctionHandleTestSuite(localfunctions);

function test_2d_grayscale %#ok<*DEFNU>

N = 2000;
img = Image.read('cameraman.tif');
sampler = RandomPositionsSampler.create(img, N);
pos = sampler.positions();

assertEqual(N, size(pos, 1));


function test_2d_color

N = 2000;
img = Image.read('peppers.png');
sampler = RandomPositionsSampler.create(img, N);
pos = sampler.positions();

assertEqual(N, size(pos, 1));


function test_3d_gray

N = 2000;
img = Image.read('brainMRI.hdr');
sampler = RandomPositionsSampler.create(img, N);
pos = sampler.positions();

assertEqual(N, size(pos, 1));

