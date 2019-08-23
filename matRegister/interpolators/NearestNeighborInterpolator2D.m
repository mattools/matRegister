classdef NearestNeighborInterpolator2D < ImageInterpolator2D
% Nearest-neighbor interpolator of a 2D image.
%   output = NearestNeighborInterpolator2D(IMG)
%
%   Example
%   I = Image2D('rice.png');
%   interp = NearestNeighborInterpolator2D(I);
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
    function obj = NearestNeighborInterpolator2D(varargin)
        % Constructs a new NearestNeighborInterpolator2D object.
        % interp = LinearInterpolator(IMG);
        % with IMG being a Image2D.
        
        if isa(varargin{1}, 'NearestNeighborInterpolator2D')
            % copy constructor
            var = varargin{1};
            image = var.Image;
        elseif isa(varargin{1}, 'Image')
            % initialisation constructor
            image = varargin{1};
            if ndims(image)~=2 %#ok<ISMAT>
                error('Image dimension should equal 2');
            end
            
        else
            error('Wrong parameter when constructing a nearest neighbor interpolator');
        end
        
        % call superclass constructor
        obj = obj@ImageInterpolator2D(image);
        
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
        
        % size of elements: number of channels by number of frames
        elSize = elementSize(obj.Image);
        
        % number of dimension of input coordinates
        dim0 = dim;
        if dim0(2) == 1
            dim0 = dim0(1);
        end
        nd = length(dim);
        
        % Create default result image
        dim2 = [dim0 elSize];
        val = ones(dim2) * obj.FillValue;
        
        % extract x and y
        xt = round(coord(:, 1));
        yt = round(coord(:, 2));
        
        % select points located inside interpolation area
        % (smaller than image physical size)
        siz = size(obj.Image);
        isInside = ~(xt < 1 | yt < 1 | xt > siz(1) | yt > siz(2));
        xt = xt(isInside);
        yt = yt(isInside);
        
        % values of the nearest neighbor
        if prod(elSize) == 1
            % case of scalar image, no movie
            val(isInside) = double(getPixels(obj.Image, xt, yt));
            
        else
            % compute interpolated values
            res = double(getPixels(obj.Image, xt, yt));
            
            % compute spatial index of each inerpolated point
            subs = cell(1, nd);
            [subs{:}] = ind2sub(dim, find(isInside));
            
            % pre-compute some indices of interpolated values
            subs2 = cell(1, length(dim2));
            subs2{nd+1} = 1:elSize(1);
            subs2{nd+2} = 1:elSize(2);
            
            % iterate on interpolated values
            for i = 1:length(subs{1})
                % compute indices of interpolated value
                for d = 1:nd
                    subInds = subs{d};
                    subs2{d} = subInds(i);
                end
                
                val(subs2{:}) = res(i, :);
            end
        end
        
        isInside = reshape(isInside, dim);
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
        isInside = ~(xt<-.5 | yt<-.5 | xt>=siz(1)-.5 | yt>=siz(2)-.5);
        xt = xt(isInside);
        yt = yt(isInside);
        isInside = reshape(isInside, dim);
        
        % indices of pixels before and after in each direction
        i1 = round(xt);
        j1 = round(yt);
        
        % values of the nearest neighbor
        val(isInside) = double(getPixels(obj.Image, i1, j1));
    end
end

end % classdef