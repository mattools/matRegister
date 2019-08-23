classdef NearestNeighborInterpolator3D < ImageInterpolator3D
% Nearest-neighbor interpolator of a 3D image.
%
%   INTERP = NearestNeighborInterpolator3D(IMG)
%
%   Example
%   I = Image2D('rice.png');
%   interp = NearestNeighborInterpolator3D(I);
%   val1 = evaluate(interp, [9.7 9.7]);
%   val2 = evaluate(interp, [10.3 10.3]);
%
%   See also
%     NearestNeighborInterpolator

% ------
% Author: David Legland
% e-mail: david.legland@inra.fr
% Created: 2010-01-07,    using Matlab 7.9.0.529 (R2009b)
% Copyright 2010 INRA - Cepia Software Platform.


%% Constructors
methods
    function obj = NearestNeighborInterpolator3D(varargin)
        % Constructs a new NearestNeighborInterpolator3D object.
        % interp = LinearInterpolator(IMG);
        % with IMG being a Image2D.
        
        if isa(varargin{1}, 'NearestNeighborInterpolator3D')
            % copy constructor
            var = varargin{1};
            image = var.Image;
        elseif isa(varargin{1}, 'Image3D')
            % initialisation constructor
            image = varargin{1};
        else
            error('Wrong parameter when constructing a nearest neighbor interpolator');
        end
        
        % call superclass constructor
        obj = obj@ImageInterpolator3D(image);
        
    end % constructor declaration
    
end % methods

methods
    function [val, isInside] = evaluate(obj, varargin)
        % Evaluate intensity of image at a given physical position.
        %
        % VAL = evaluate(INTERP, POS);
        % where POS is a Nx2 array containing alues of x- and y-coordinates
        % of positions to evaluate image, return an array with as many
        % values as POS.
        %
        % VAL = evaluate(INTERP, X, Y)
        % X and Y should be the same size. The result VAL has the same size
        % as X and Y.
        %
        % [VAL, INSIDE] = evaluate(...)
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
        
        % number of positions to process
        N = size(coord, 1);
        
        % extract x and y
        xt = coord(:, 1);
        yt = coord(:, 2);
        zt = coord(:, 3);
        
        % select points located inside interpolation area
        % (smaller than image physical size)
        siz = size(obj.Image);
        isBefore    = sum(coord <  .5, 2) > 0;
        isAfter     = sum(coord >= (siz(ones(N,1), :))+.5, 2) > 0;
        isInside    = ~(isBefore | isAfter);
        
        xt = xt(isInside);
        yt = yt(isInside);
        zt = zt(isInside);
        isInside = reshape(isInside, dim);
        
        % indices of pixels before and after in each direction
        i1 = round(xt);
        j1 = round(yt);
        k1 = round(zt);
        
        % values of the nearest neighbor
        val(isInside) = double(getPixels(obj.Image, i1, j1, k1));
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
        
        % indices of pixels before and after in each direction
        i1 = round(xt);
        j1 = round(yt);
        k1 = round(zt);
        
        % values of the nearest neighbor
        val(isInside) = double(getPixels(obj.Image, i1, j1, k1));
    end
end

end % classdef