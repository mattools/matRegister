function fval = computeValue(obj)
%COMPUTEVALUE Compute metric value using current state
%
%   FVAL = obj.computeValue()
%
%   Example
%   computeValue
%
%   See also
%

% ------
% Author: David Legland
% e-mail: david.legland@inra.fr
% Created: 2011-01-06,    using Matlab 7.9.0.529 (R2009b)
% Copyright 2011 INRA - Cepia Software Platform.


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
    
    res(i) = computeMeanSquaredDifferences(...
        obj.Images{i1}, obj.Images{i2}, obj.Points);
end

% sum of SSD computed over couples
fval = sum(res);


function res = computeMeanSquaredDifferences(img1, img2, points)

% compute values in image 1
values1 = evaluate(img1, points);

% compute values in image 2
values2 = evaluate(img2, points);


% compute squared differences
diff = (values2 - values1).^2;

% Sum of squared differences normalized by number of test points
res = mean(diff);
