%TESTSAMPLERS  Basic test cases for image samplers
%
%   output = testSamplers(input)
%
%   Example
%   testSamplers
%
%   See also
%
%
% ------
% Author: David Legland
% e-mail: david.legland@grignon.inra.fr
% Created: 2011-07-20,    using Matlab 7.9.0.529 (R2009b)
% Copyright 2011 INRA - Cepia Software Platform.



img = Image.read('cameraman.tif');

sampler = FullImageSampler.create(img);

pos = sampler.positions();

assertEqual(prod(size(img)), size(pos, 1)); %#ok<PSIZE>
