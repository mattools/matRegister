classdef LinearInterpolator3D < ImageInterpolator3D
%LINEARINTERPOLATOR3D  Linear interpolator of an image.
%   INTERP = LinearInterpolator3D(IMG)
%
%   Example
%   
%   See also
%

% ------
% Author: David Legland
% e-mail: david.legland@inra.fr
% Created: 2010-01-07,    using Matlab 7.9.0.529 (R2009b)
% Copyright 2010 INRA - Cepia Software Platform.


%% Constructors
methods
    function obj = LinearInterpolator3D(varargin)
        % Constructs a new LinearInterpolator object.
        % interp = LinearInterpolator(IMG);
        % with IMG being a Image2D.
        
        % process input arguments
        if isa(varargin{1}, 'LinearInterpolator3D')
            % copy constructor
            var = varargin{1};
            img = var.Image;
        elseif isa(varargin{1}, 'Image')
            % copy constructor
            img = varargin{1};
        else
            error('Wrong parameter when constructing a 3D linear interpolator');
        end

        % call superclass constructor
        obj = obj@ImageInterpolator3D(img);

    end % constructor declaration    
end % methods

methods
    function d = getDimension(obj) %#ok<MANU>
        % Dimension of the interpolated image.
        d = 3;
    end
    
    function [value, isInside] = evaluate(obj, varargin)
        % Evaluate intensity of image at a given physical position.
        %
        % VAL = evaluate(INTERP, POS);
        % where POS is a N-by-2 array containing alues of x- and y-coordinates
        % of positions to evaluate in image, return an array with as many values as
        % POS.
        %
        % VAL = evaluate(INTERP, X, Y)
        % where X and Y are two arrays the same size, return the result in an array
        % VAL that has the same size as X and Y.
        %
        % [VAL, INSIDE] = evaluate(...)
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
    end
    
    function [val, isInside] = evaluateAtIndex(obj, varargin)
        % Evaluate the value of an image for a point given in image coord.
        
        % eventually convert inputs to a single nPoints-by-ndims array
        [index, dim] = ImageFunction.mergeCoordinates(varargin{:});
        
        % number of positions to process
        N = size(index, 1);
        
        % initialize result with default value
        val = ones(dim) * obj.FillValue;
        
        % extract x and y
        xt = index(:, 1);
        yt = index(:, 2);
        zt = index(:, 3);
        
        % select points located inside interpolation area
        % (smaller than image size)
        siz = size(obj.Image);
        isBefore    = sum(index<1, 2)<0;
        isAfter     = sum(index>=(siz(ones(N,1), :)), 2)>0;
        isInside = ~(isBefore | isAfter);
        isInside = reshape(isInside, dim);
        
        % keep only valid indices for computation
        xt = xt(isInside);
        yt = yt(isInside);
        zt = zt(isInside);
        
        % compute distances to lower-left pixel
        i1 = floor(xt);
        j1 = floor(yt);
        k1 = floor(zt);
        
        % calcule les distances au bord inferieur bas et gauche
        dx = (xt-i1);
        dy = (yt-j1);
        dz = (zt-k1);
        
        % values of the 4 pixels around each point
        val111 = double(getPixels(obj.Image, i1, j1, k1)).*(1-dx).*(1-dy).*(1-dz);
        val211 = double(getPixels(obj.Image, i1, j1+1, k1)).*(1-dx).*dy.*(1-dz);
        val121 = double(getPixels(obj.Image, i1+1, j1, k1)).*dx.*(1-dy).*(1-dz);
        val221 = double(getPixels(obj.Image, i1+1, j1+1, k1)).*dx.*dy.*(1-dz);
        val112 = double(getPixels(obj.Image, i1, j1, k1+1)).*(1-dx).*(1-dy).*dz;
        val212 = double(getPixels(obj.Image, i1, j1+1, k1+1)).*(1-dx).*dy.*dz;
        val122 = double(getPixels(obj.Image, i1+1, j1, k1+1)).*dx.*(1-dy).*dz;
        val222 = double(getPixels(obj.Image, i1+1, j1+1, k1+1)).*dx.*dy.*dz;
        
        % compute result values
        val(isInside) = ...
            val111 + val121 + val211 + val221 + ...
            val112 + val122 + val212 + val222;
    end
end

end % classdef