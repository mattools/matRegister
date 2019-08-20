function [value, isInside] = evaluate(obj, varargin)
% Evaluate intensity of image at a given physical position.
%
% VAL = INTERP.evaluate(POS);
% where POS is a N-by-2 array containing alues of x- and y-coordinates
% of positions to evaluate in image, return an array with as many values as
% POS. 
%
% VAL = INTERP.evaluate(X, Y)
% where X and Y are two arrays the same size, return the result in an array
% VAL that has the same size as X and Y. 
%
% [VAL INSIDE] = INTERP.evaluate(...)
% Also return a boolean flag the same size as VAL indicating whether or not
% the given position as located inside the evaluation frame.
%


% eventually convert inputs to a single nPoints-by-ndims array
[point, dim] = ImageFunction.mergeCoordinates(varargin{:});

% Evaluates image value for a given position
coord = pointToContinuousIndex(obj.Image, point);

% extract x and y
xt  = coord(:, 1);
yt  = coord(:, 2);
zt  = coord(:, 3);

% select points located inside interpolation area
% (smaller than image physical size)
siz = size(obj.Image);
isBefore    = sum(coord < 1, 2) > 0;
isAfter     = xt >= siz(1) | yt >= siz(2) | zt >= siz(3);
isInside    = ~(isBefore | isAfter);
isInside    = reshape(isInside, dim);

% keep only valid indices for computation
xt  = xt(isInside);
yt  = yt(isInside);
zt  = zt(isInside);

% indices of pixels before in each direction
i1  = floor(xt);
j1  = floor(yt);
k1  = floor(zt);

% compute distances to lower-left pixel
dx  = xt - i1;
dy  = yt - j1;
dz  = zt - k1;

dxi = 1 - dx;
dyi = 1 - dy;
dzi = 1 - dz;

% image sizes
dimX    = siz(1);
dimXY   = siz(1) * siz(2);

% pre-compute indices for linear indexing
inds    = (k1 - 1) * dimXY + (j1 - 1) * dimX + i1;

% values of the 4 pixels around each point
vals = double(obj.Image.Data(inds))                 .*dxi   .* dyi  .* dzi;
vals = vals + double(obj.Image.Data(inds + 1))      .*dx    .* dyi  .* dzi;
vals = vals + double(obj.Image.Data(inds + dimX))   .*dxi   .* dy   .* dzi;
vals = vals + double(obj.Image.Data(inds + dimX+1)) .*dx    .* dy   .* dzi;
inds = inds + dimXY;
vals = vals + double(obj.Image.Data(inds))          .* dxi  .* dyi  .* dz;
vals = vals + double(obj.Image.Data(inds + 1))      .* dx   .* dyi  .* dz;
vals = vals + double(obj.Image.Data(inds + dimX))   .* dxi  .* dy   .* dz;
vals = vals + double(obj.Image.Data(inds + dimX+1)) .* dx   .* dy   .* dz;

% Create default result image
value = zeros(dim);
value(:) = obj.FillValue;
value(isInside) = vals;
