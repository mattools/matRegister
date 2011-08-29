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
%
% ------
% Author: David Legland
% e-mail: david.legland@grignon.inra.fr
% Created: 2011-03-15,    using Matlab 7.9.0.529 (R2009b)
% Copyright 2011 INRA - Cepia Software Platform.


%% Declaration of class properties
properties
    % The size of the resulting image
    outputSize = [];

    % Data type of resulting image, default is 'uint8'.
    outputType = 'double';
    
    % physical origin of image first point
    origin = [];

    % Physical spacing between pixels
    spacing = [];
end

%TODO: specify output type;

%% Constructors
methods
    function this = TransformJacobianResampler(varargin)
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
            this.outputSize = var.outputSize;
            this.origin     = var.origin;
            this.spacing    = var.spacing;
            
        elseif isa(varargin{1}, 'Image')
            % Initialize fields from a given image
            img = varargin{1};
            this.outputSize = size(img);
            this.origin     = img.getOrigin();
            this.spacing    = img.getSpacing();
            
        elseif isnumeric(varargin{1})
            % Gives lx and ly
            nbDims = length(varargin);
            this.outputSize = zeros(1, nbDims);
            this.origin     = zeros(1, nbDims);
            this.spacing    = ones(1, nbDims);
            for i=1:nbDims
                var = varargin{i};
                this.outputSize(i)  = length(var);
                this.origin(i)      = var(1);
                this.spacing(i)     = var(2)-var(1);
            end
        else
            error('Wrong parameter when constructing an ImageResampler');
        end
        
    end % constructor declaration
end % methods

methods
    function setOutputType(this, dataType)
        % Specifiy the type of resulting image
        this.outputType = dataType;
    end
    
    function type = getOutputType(this)
        % Return the type of resulting image
        type = this.outputType;
    end
    
    function img2 = resample(this, transform)
        % Resample a transform
        %
        
        if nargin < 2
            error('Need to specify a transform');
        end
        
        if ~isa(transform, 'Transform')
            error('First argument must be a transform');
        end
        
        % precompute grid basis for 2D
        lx = (0:this.outputSize(1)-1)*this.spacing(1) + this.origin(1);
        ly = (0:this.outputSize(2)-1)*this.spacing(2) + this.origin(2);
        
        % different processing depending on image dimension
        outputDim = length(this.outputSize);
        if outputDim == 2
            % Process 2D images
            
            % sample grid
            [x y] = meshgrid(lx, ly);

            % initialize result array
            vals = zeros(size(x), this.outputType);
        
            % compute all jacobians (result stored in a 2*2*N array)
            jac = getJacobian(transform, [x(:) y(:)]);
            for i = 1:numel(x)
                vals(i) = det(jac(:,:,i));
            end

        elseif outputDim == 3
            % Process 3D images
            
            % sample grid
            lz = (0:this.outputSize(3)-1)*this.spacing(3) + this.origin(3);
            [x y z] = meshgrid(lx, ly, lz);
            
            % initialize result array
            vals = zeros(size(x), this.outputType);
        
            % compute result values
            try
                % try first with full call
                jac = transform.getJacobian([x(:) y(:) z(:)]);
                for i = 1:numel(x)
                    vals(i) = det(jac(:,:,i));
                end
                
            catch 
                % if memory limit is reached, performs element by element
                warning([mfilename ':MemoryLimit'], ...
                    'Memory limit reached - switched to elementwise computation');
                for i = 1:numel(x)
                    jac = transform.getJacobian([x(i) y(i) z(i)]);
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
        img2.setOrigin(this.origin);
        img2.setSpacing(this.spacing);
    end

end % methods


end % classef 
