function test_suite = test_RandomSampler(varargin)
%TEST_RANDOMSAMPLER  One-line description here, please.
%
%   output = test_RandomSampler(input)
%
%   Example
%   test_RandomSampler
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
sampler = RandomSampler.create(img, N);
pos = sampler.positions();

assertEqual(N, size(pos, 1));


function test_2d_color

N = 2000;
img = Image.read('peppers.png');
sampler = RandomSampler.create(img, N);
pos = sampler.positions();

assertEqual(N, size(pos, 1));


function test_3d_gray

N = 2000;
img = Image.read('brainMRI.hdr');
sampler = RandomSampler.create(img, N);
pos = sampler.positions();

assertEqual(N, size(pos, 1));

