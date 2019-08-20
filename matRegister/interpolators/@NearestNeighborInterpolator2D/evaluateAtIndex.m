function [val, isInside] = evaluateAtIndex(obj, varargin)
% Evaluate the value of an image for a point given in image coord.

% eventually convert inputs to a single nPoints-by-ndims array
[index, dim] = ImageFunction.mergeCoordinates(varargin{:});

val = ones(dim) * obj.FillValue;

% extract x and y
xt = index(:, 1);
yt = index(:, 2);

% select points located inside interpolation area
% (smaller than image size)
siz = size(obj.Image);
isInside = ~(xt<-.5 | yt<-.5 | xt>=siz(1)-.5 | yt>=siz(2)-.5);
xt = xt(isInside);
yt = yt(isInside);
isInside = reshape(isInside, dim);

% indices of pixels before and after in each direction
i1 = round(xt);
j1 = round(yt);

% values of the nearest neighbor
val(isInside) = double(getPixels(obj.Image, i1, j1));
