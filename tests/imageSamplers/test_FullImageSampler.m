function test_suite = test_FullImageSampler(varargin) %#ok<STOUT>
%TEST_FULLIMAGESAMPLER  One-line description here, please.
%
%   output = test_FullImageSampler(input)
%
%   Example
%   test_FullImageSampler
%
%   See also
%
%
% ------
% Author: David Legland
% e-mail: david.legland@grignon.inra.fr
% Created: 2011-07-20,    using Matlab 7.9.0.529 (R2009b)
% Copyright 2011 INRA - Cepia Software Platform.

initTestSuite;

function test_2d_grayscale %#ok<*DEFNU>

img = Image.read('cameraman.tif');
sampler = FullImageSampler.create(img);
pos = sampler.positions();

assertEqual(prod(size(img)), size(pos, 1)); %#ok<PSIZE>


function test_2d_color

img = Image.read('peppers.png');
sampler = FullImageSampler.create(img);
pos = sampler.positions();

assertEqual(prod(size(img)), size(pos, 1)); %#ok<PSIZE>


function test_3d_gray

img = Image.read('brainMRI.hdr');
sampler = FullImageSampler.create(img);
pos = sampler.positions();

assertEqual(prod(size(img)), size(pos, 1)); %#ok<PSIZE>

