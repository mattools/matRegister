%TESTIMAGE2D  One-line description here, please.
%   output = testImage2D(input)
%
%   Example
%   testImage2D
%
%   See also
%
%
% ------
% Author: David Legland
% e-mail: david.legland@grignon.inra.fr
% Created: 2009-04-20,    using Matlab 7.7.0.471 (R2008b)
% Copyright 2009 INRA - Cepia Software Platform.
% Licensed under the terms of the LGPL, see the file "license.txt"

% clean up
clear all;
clear classes;

% create test image with 4 columns and 3 rows
dat = [1 2 3 4;5 6 7 8;9 10 11 12];
img = Image.create(dat);

% get a pixel
assertEqual(11, img(3, 3));

% get a pixel not on diagonal
assertEqual(9, img(1, 3));

% get size
assertElementsAlmostEqual([4 3], size(img));

