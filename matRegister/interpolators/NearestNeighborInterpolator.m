classdef NearestNeighborInterpolator < ImageInterpolator
% Nearest-neighbor interpolator of an image.
%
%   INTERP = NearestNeighborInterpolator(IMG)
%
%   Example
%   I = Image.read('rice.png');
%   interp = NearestNeighborInterpolator(I);
%   val1 = interp.evaluate([9.7 9.7]);
%   val2 = interp.evaluate([10.3 10.3]);
%
%   See also
%     ImageInterpolator
%

% ------
% Author: David Legland
% e-mail: david.legland@inra.fr
% Created: 2010-01-07,    using Matlab 7.9.0.529 (R2009b)
% Copyright 2010 INRA - Cepia Software Platform.


%% Constructors
methods
    function obj = NearestNeighborInterpolator(varargin)
        % Constructs a new NearestNeighborInterpolator object.
        % interp = LinearInterpolator(IMG);
        % with IMG being a Image2D.
        
        if isa(varargin{1}, 'NearestNeighborInterpolator')
            % copy constructor
            var = varargin{1};
            image = var.Image;
        elseif isa(varargin{1}, 'Image')
            % initialisation constructor
            image = varargin{1};
                        
        else
            error('Wrong parameter when constructing a nearest neighbor interpolator');
        end
        
        % call superclass constructor
        obj = obj@ImageInterpolator(image);
        
    end % constructor declaration
    
end % methods

%% Methods
methods
    function [val, isInside] = evaluate(obj, varargin)
        % Evaluate intensity of image at a given physical position.
        %
        %   VAL = evaluate(INTERP, POS);
        %   where POS is a N-by-2 or N-by-3 array containing values of x-, y-, and
        %   eventually z-coordinates of positions to evaluate image, returns a
        %   column array with as many rows as the number of rows in POS.
        %
        %   VAL = evaluate(INTERP, X, Y)
        %   VAL = evaluate(INTERP, X, Y, Z)
        %   X, Y and Z should be the same size. The result VAL has the same size as
        %   X and Y.
        %
        %
        %   [VAL INSIDE] = evaluate(INTERP, ...)
        %   Also return a boolean flag the same size as VAL indicating whether or
        %   not the given position as located inside the evaluation frame.
        %
        
        % number of dimensions of base image
        nd = ndims(obj.Image);
        
        % size of elements: number of channels by number of frames
        elSize = elementSize(obj.Image);
        
        
        % eventually convert inputs to a single nPoints-by-ndims array
        [point, dim] = ImageFunction.mergeCoordinates(varargin{:});
        
        if size(point, 2) ~= nd
            error('Dimension of input positions should be the same as image');
        end
        
        % Evaluates image value for a given position
        coord = pointToContinuousIndex(obj.Image, point);
        
        
        % size and number of dimension of input coordinates
        dim0 = dim;
        if dim0(end) == 1
            dim0 = dim0(1:end-1);
        end
        resNDim = length(dim0);
        
        % Create default result image
        val = zeros([dim0 elSize]);
        val(:) = obj.FillValue;
        
        % extract x and y
        xt = coord(:, 1);
        yt = coord(:, 2);
        if nd > 2
            zt = coord(:, 3);
        end
        
        % number of positions to process
        N = size(coord, 1);
        
        % select points located inside interpolation area
        % (smaller than image physical size)
        siz = size(obj.Image);
        isBefore    = sum(coord < .5, 2)>0;
        isAfter     = sum(coord >= (siz(ones(N,1), :))+.5, 2)>0;
        isInside    = ~(isBefore | isAfter);
        
        xt = xt(isInside);
        yt = yt(isInside);
        isInside = reshape(isInside, dim);
        
        % indices of pixels before and after in each direction
        i1 = round(xt);
        j1 = round(yt);
        if nd > 2
            zt = zt(isInside);
            k1 = round(zt);
        end
        
        % values of the nearest neighbor
        if prod(elSize) == 1
            % case of scalar image elements (no vector image, no movie image)
            if nd == 2
                val(isInside) = double(getPixels(obj.Image, i1, j1));
            else
                val(isInside) = double(getPixels(obj.Image, i1, j1, k1));
            end
        else
            % If dimension number of elements is >1, need more processing.
            
            % compute interpolated values
            if nd == 2
                res = double(getPixels(obj.Image, i1, j1));
            else
                res = double(getPixels(obj.Image, i1, j1, k1));
            end
            
            % compute spatial index of each interpolated point
            subs = cell(1, resNDim);
            [subs{:}] = ind2sub(dim, find(isInside));
            
            % pre-compute some indices of interpolated values
            subs2 = num2cell(ones(1, length(elSize)));
            subs2{resNDim+1} = 1:elSize(1);
            subs2{resNDim+2} = 1:elSize(2);
            
            % iterate on interpolated values
            for i = 1:length(subs{1})
                % compute matrix indices of interpolated value
                for d = 1:resNDim
                    subInds = subs{d};
                    subs2{d} = subInds(i);
                end
                
                % places components of interpolated value at the right loaction
                val(subs2{:}) = res(i, :);
            end
        end
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