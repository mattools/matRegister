function val = computeValue(obj)
%COMPUTEVALUE Compute metric value 
%
%   output = computeValue(input)
%
%   Example
%   computeValue
%
%   See also
%

% ------
% Author: David Legland
% e-mail: david.legland@inra.fr
% Created: 2010-09-29,    using Matlab 7.9.0.529 (R2009b)
% Copyright 2010 INRA - Cepia Software Platform.


% extract number of images
nImages = length(obj.Images);

% generate all possible couples of images
combis  = sortrows(combnk(1:nImages, 2));
nCombis = size(combis, 1);

% compute SSD for each image couple
res = zeros(nCombis, 1);
for i = 1:nCombis
    i1 = combis(i,1);
    i2 = combis(i,2);
    
    res(i) = computeSSDMetric(obj.Images{i1}, obj.Images{i2}, obj.Points);
end

% sum of SSD computed over couples
val = sum(res);

function res = computeSSDMetric(img1, img2, points)

% compute values in image 1
[values1, inside1] = evaluate(img1, points);

% compute values in image 2
[values2, inside2] = evaluate(img2, points);

% keep only valid values
inds = inside1 & inside2;

% compute result
diff = (values2(inds) - values1(inds)).^2;
res = sum(diff);

