classdef ImageResampler < handle
% Resample an image using a given spatial basis.
%
%   RES = ImageResampler(LX, LY);
%   RES = ImageResampler(LX, LY, LZ);
%   RES = ImageResampler(BASE_IMAGE);
%
% %   the following is not yet implemented
% %   RES = ImageResampler('outputSize', [100 100], 'origin', [0 0], ...
% %       'spacing', [1 1]);
%
%   Usage:
%   IMG = Image2D('cameraman.tif');
%   IMG2 = RES.resample(IMG, 'linear');
%
%   Example
%   ImageResampler
%
%   See also
%

% ------
% Author: David Legland
% e-mail: david.legland@inra.fr
% Created: 2010-06-03,    using Matlab 7.9.0.529 (R2009b)
% Copyright 2010 INRA - Cepia Software Platform.

%% Declaration of class properties
properties
    % The size of the resulting image
    OutputSize = [];

    % Data type of resulting image, default is 'uint8'.
    OutputType = 'uint8';
    
    % physical origin of image first point
    Origin = [];

    % Physical spacing between pixels
    Spacing = [];
end

%TODO: specify output type;

%% Constructors
methods
    function obj = ImageResampler(varargin)
        % Constructs a new ImageResampler object.

        if isa(varargin{1}, 'ImageResampler')
            % copy constructor

            % copy each field
            var = varargin{1};
            obj.OutputSize = var.OutputSize;
            obj.OutputType = var.OutputType;
            obj.Origin     = var.Origin;
            obj.Spacing    = var.Spacing;
            
        elseif isa(varargin{1}, 'Image')
            % Initialize fields from a given image
            img = varargin{1};
            obj.OutputSize = size(img);
            obj.OutputType = img.getDataType();
            obj.Origin     = img.Origin;
            obj.Spacing    = img.Spacing;
            
        elseif isnumeric(varargin{1})
            % Gives lx and ly
            nbDims = length(varargin);
            obj.OutputSize = zeros(1, nbDims);
            obj.Origin     = zeros(1, nbDims);
            obj.Spacing    = ones(1, nbDims);
            for i=1:nbDims
                var = varargin{i};
                obj.OutputSize(i)  = length(var);
                obj.Origin(i)      = var(1);
                obj.Spacing(i)     = var(2)-var(1);
            end
        else
            error('Wrong parameter when constructing an ImageResampler');
        end
    end % constructor declaration
end % methods

methods
    function setOutputType(obj, dataType)
        % Specifiy the type of resulting image
        obj.OutputType = dataType;
    end
    
    function type = getOutputType(obj)
        % Return the type of resulting image
        type = obj.OutputType;
    end
    
    function img2 = resample(obj, varargin)
        % Resample an image, or an interpolated image
        %
        % IMG2 = resample(OBJ, IMGFUN);
        % Use image function IMGFUN as reference. IMGFUN is an instance of 
        % ImageFunction, such as ImageInterpolator.
        %
        % IMG2 = resample(OBJ, IMG);
        % Resample image IMG using linear interpolation.
        %
        % IMG2 = resample(OBJ, IMG, TYPE);  
        % Resample image with specified interpolation type. TYPE can be:
        % 'linear' (default)
        % 'nearest', but is only supported for 2D images.
        %
        
        if isempty(varargin)
            error('Need to specify an image to resample');
        end
        
        var = varargin{1};
        if isa(var, 'ImageFunction')
            % if input is already an ImageFunction, nothing to do
            imgFun = var;
        else
            % extract input image
            if isa(var, 'Image')
                img = var;
            elseif isnumeric(var)
                img = Image.create(var);
            else
                error('Please specify image to resample');
            end
            
            % determines interpolation type
            interpolationType = 'linear';
            if length(varargin)>1
                interpolationType = varargin{2};
            end
            
            % create the appropriate image interpolator
            imgFun = ImageInterpolator.create(img, interpolationType);
        end
        
        
        % precompute positions for 2D
        lx = (0:obj.OutputSize(1)-1)*obj.Spacing(1) + obj.Origin(1);
        ly = (0:obj.OutputSize(2)-1)*obj.Spacing(2) + obj.Origin(2);
        
        % different processing depending on image dimension
        outputDim = length(obj.OutputSize);
        if outputDim==2
            % Process 2D images
            
            [x, y] = meshgrid(lx, ly);

            vals = zeros(size(x), obj.OutputType);
        
            vals(:) = evaluate(imgFun, [x(:) y(:)]);
            
        elseif outputDim==3
            % Process 3D images
            
            lz = (0:obj.OutputSize(3)-1)*obj.Spacing(3) + obj.Origin(3);
            [x, y, z] = meshgrid(lx, ly, lz);
            
            vals = zeros(size(x), obj.OutputType);
        
            vals(:) = evaluate(imgFun, [x(:) y(:) z(:)]);
            
        else
            % TODO: implement for greater dimensions ?
            error('Resampling is only implemented for dimensions 2 and 3');
        end
        
        % create new image from computed values
        img2 = Image.create(vals);
        
        % copy spatial calibration info to image
        img2.Origin = obj.Origin;
        img2.Spacing = obj.Spacing;
    end

end % abstract methods

end  % classdef
