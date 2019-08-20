classdef TransformJacobianResampler < handle
%TRANSFORMJACOBIANRESAMPLER Create a new image showing jacobian of a transform 
%
%   Resample the jacobian matrix determinant of a transform, and build a
%   grayscale image to represent it.
%
%   Example
%     TJR = TransformJacobianResampler(REFIMG);
%     JACIMG = resample(TJR, TRANSFO);
%     show(JACIMG)
%
%   See also
%     ImageResampler
%

% ------
% Author: David Legland
% e-mail: david.legland@inra.fr
% Created: 2011-03-15,    using Matlab 7.9.0.529 (R2009b)
% Copyright 2011 INRA - Cepia Software Platform.


%% Declaration of class properties
properties
    % The size of the resulting image
    OutputSize = [];

    % Data type of resulting image, default is 'uint8'.
    OutputType = 'double';
    
    % physical origin of image first point
    Origin = [];

    % Physical spacing between pixels
    Spacing = [];
end

%TODO: specify output type;

%% Constructors
methods
    function obj = TransformJacobianResampler(varargin)
        % Constructs a new TransformJacobianResampler object.
        %
        % TJR = TransformJacobianResampler(REFIMG)
        % Creates a resampler using the same spatial basis as the given
        % reference image.
        %
        % TJR = TransformJacobianResampler(TJR0)
        % Copy constructor
        %
        %

        if isempty(varargin)
            error('TransformJacobianResampler:WrongNumberOfARguments', ...
                'Need to specify at least one argument to initialise TransformJacobianResampler');
        end
        
        if isa(varargin{1}, 'TransformJacobianResampler')
            % copy constructor

            % copy each field
            var = varargin{1};
            obj.OutputSize = var.OutputSize;
            obj.Origin     = var.Origin;
            obj.Spacing    = var.Spacing;
            
        elseif isa(varargin{1}, 'Image')
            % Initialize fields from a given image
            img = varargin{1};
            obj.OutputSize = size(img);
            obj.Origin     = img.getOrigin();
            obj.Spacing    = img.getSpacing();
            
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
    
    function img2 = resample(obj, transform)
        % Resample a transform
        %
        
        if nargin < 2
            error('Need to specify a transform');
        end
        
        if ~isa(transform, 'Transform')
            error('First argument must be a transform');
        end
        
        % precompute grid basis for 2D
        lx = (0:obj.OutputSize(1)-1)*obj.Spacing(1) + obj.Origin(1);
        ly = (0:obj.OutputSize(2)-1)*obj.Spacing(2) + obj.Origin(2);
        
        % different processing depending on image dimension
        outputDim = length(obj.OutputSize);
        if outputDim == 2
            % Process 2D images
            
            % sample grid
            [x, y] = meshgrid(lx, ly);

            % initialize result array
            vals = zeros(size(x), obj.OutputType);
        
            % compute all jacobians (result stored in a 2*2*N array)
            jac = jacobianMatrix(transform, [x(:) y(:)]);
            for i = 1:numel(x)
                vals(i) = det(jac(:,:,i));
            end

        elseif outputDim == 3
            % Process 3D images
            
            % sample grid
            lz = (0:obj.OutputSize(3)-1)*obj.Spacing(3) + obj.Origin(3);
            [x, y, z] = meshgrid(lx, ly, lz);
            
            % initialize result array
            vals = zeros(size(x), obj.OutputType);
        
            % compute result values
            try
                % try first with full call
                jac = jacobianMatrix(transform, [x(:) y(:) z(:)]);
                for i = 1:numel(x)
                    vals(i) = det(jac(:,:,i));
                end
                
            catch 
                % if memory limit is reached, performs element by element
                warning([mfilename ':MemoryLimit'], ...
                    'Memory limit reached - switched to element-wise computation');
                for i = 1:numel(x)
                    jac = jacobianMatrix(transform, [x(i) y(i) z(i)]);
                    vals(i) = det(jac);
                end
            end
            
        else
            % TODO: implement for greater dimensions ?
            error('Resampling is only implemented for dimensions 2 and 3');
        end
        
        % create a new Image from the result
        img2 = Image.create(vals);
        
        % copy spatial calibration info to image
        img2.Origin = obj.Origin;
        img2.Spacing = obj.Spacing;
    end

end % methods


end % classef 
