classdef LinearInterpolator2D < ImageInterpolator2D
%LINEARINTERPOLATOR2D  Linear interpolator of a 2D image.
%   output = LinearInterpolator2D(input)
%
%   Example
%   I = Image2D('rice.png');
%   interp = LinearInterpolator2D(I);
%   val00 = I.getPixel(10, 10);
%   val01 = I.getPixel(11, 10);
%   val10 = I.getPixel(10, 11);
%   val11 = I.getPixel(11, 11);
%   valInt = interp.evaluate([10.5 10.5]);
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
    function obj = LinearInterpolator2D(varargin)
        % Constructs a new LinearInterpolator2D object.
        % interp = LinearInterpolator(IMG);
        % with IMG being a Image2D.
        
        if isa(varargin{1}, 'LinearInterpolator2D')
            % copy constructor
            var = varargin{1};
            img = var.Image;
        elseif isa(varargin{1}, 'Image')
            % copy constructor
            img = varargin{1};
        else
            error('Wrong parameter when constructing a linear interpolator');
        end
        
        obj = obj@ImageInterpolator2D(img);
    end % constructor declaration    
end % methods

methods
    function d = getDimension(obj) %#ok<MANU>
        % Dimension of the interpolated image.
        d = 2;
    end

    function [val, isInside] = evaluate(obj, varargin)
        % Evaluate intensity of image at a given physical position.
        %
        % VAL = INTERP.evaluate(POS);
        % where POS is a Nx2 array containing values of x- and y-coordinates
        % of positions to evaluate image, return an array with as many
        % values as POS.
        %
        % VAL = INTERP.evaluate(X, Y)
        % X and Y should be the same size. The result VAL has the same size
        % as X and Y.
        %
        % [VAL INSIDE] = INTERP.evaluate(...)
        % Also return a boolean flag the same size as VAL indicating
        % whether or not the given position as located inside the
        % evaluation frame.
        %
        
        
        % eventually convert inputs to a single nPoints-by-ndims array
        [point, dim] = ImageFunction.mergeCoordinates(varargin{:});
        
        % Evaluates image value for a given position
        coord = pointToContinuousIndex(obj.Image, point);
        
        % Create default result image
        val = ones(dim) * obj.FillValue;
        
        % extract x and y
        xt  = coord(:, 1);
        yt  = coord(:, 2);
        
        % select points located inside interpolation area
        % (smaller than image physical size)
        siz = size(obj.Image);
        isInside = ~(xt < 1 | yt < 1 | xt >= siz(1) | yt >= siz(2));
        xt  = xt(isInside);
        yt  = yt(isInside);
        isInside = reshape(isInside, dim);
        
        % indices of pixels before and after in each direction
        i1  = floor(xt);
        j1  = floor(yt);
        
        % compute distances to lower-left pixel
        dx  = xt - i1;
        dy  = yt - j1;
        
        % pre-compute distance to next pixel
        dxi = 1 - dx;
        dyi = 1 - dy;
        
        % image sizes
        siz     = size(obj.Image);
        dimX    = siz(1);
        
        inds    = (j1 - 1) * dimX + i1;
        
        % values of the 4 pixels around each point
        val11 = double(obj.Image.Data(inds))          .*dxi   .* dyi;
        val12 = double(obj.Image.Data(inds + 1))      .*dx    .* dyi;
        val21 = double(obj.Image.Data(inds + dimX))   .*dxi   .* dy;
        val22 = double(obj.Image.Data(inds + dimX+1)) .*dx    .* dy;

        % compute result values
        val(isInside) = val11 + val12 + val21 + val22;
    end
    
    
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
        isInside = ~(xt<1 | yt<1 | xt>=siz(1) | yt>=siz(2));
        xt = xt(isInside);
        yt = yt(isInside);
        isInside = reshape(isInside, dim);
        
        % compute distances to lower-left pixel
        i1 = floor(xt);
        j1 = floor(yt);
        
        % calcule les distances au bord inferieur bas et gauche
        dx = (xt-i1);
        dy = (yt-j1);
        
        % values of the 4 pixels around each point
        val11 = double(getPixels(obj.Image, i1, j1)).*(1-dx).*(1-dy);
        val21 = double(getPixels(obj.Image, i1, j1+1)).*(1-dx).*dy;
        val12 = double(getPixels(obj.Image, i1+1, j1)).*dx.*(1-dy);
        val22 = double(getPixels(obj.Image, i1+1, j1+1)).*dx.*dy;
        
        % compute result values
        val(isInside) = val11 + val12 + val21 + val22;
    end
end

end % classdef
