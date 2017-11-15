function test_suite = test_GridImageSampler(varargin) 
%TEST_GRIDIMAGESAMPLER  One-line description here, please.
%
%   output = test_GridImageSampler(input)
%
%   Example
%   test_GridImageSampler
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

img = Image.read('cameraman.tif');
steps = [8 8];
sampler = GridImageSampler.create(img, steps);
pos = sampler.positions();

gridSize = floor(size(img) ./ steps);
assertEqual(prod(gridSize), size(pos, 1));


function test_2d_color

img = Image.read('peppers.png');
steps = [8 8];
sampler = GridImageSampler.create(img, steps);
pos = sampler.positions();

gridSize = floor(size(img) ./ steps);
assertEqual(prod(gridSize), size(pos, 1));


function test_3d_gray

img = Image.read('brainMRI.hdr');
steps = [8 8 2];
sampler = GridImageSampler.create(img, steps);
pos = sampler.positions();

gridSize = floor((size(img) - 1) ./ steps) + 1;
assertEqual(prod(gridSize), size(pos, 1));

