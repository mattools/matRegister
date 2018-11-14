classdef CenteredMotionTransform2D < AffineTransform & ParametricTransform & CenteredTransformAbstract
%Transformation model for a centered rotation followed by a translation
%   
%   Inner optimisable parameters of the transform have the following form:
%   params[1] = theta, in degrees
%   params[2] = tx
%   params[3] = ty
% 
%   TRANS = CenteredMotionTransform2D();
%   Create using default parameters (all zero).
%
%   TRANS = CenteredMotionTransform2D(PARAMS);
%   Create a new transform model by initializing parameters. PARAMS is a
%   1-by-3 row vector containing TX, TY, and THETA parameters.
%
%   TRANS = CenteredMotionTransform2D('center', CENTER);
%   Initialize the center of the transform. CENTER must be a 1-by-2 row
%   vector containing coordinates of the transform center.
%

% ------
% Author: David Legland
% e-mail: david.legland@inra.fr
% Created: 2010-02-17,    using Matlab 7.9.0.529 (R2009b)
% Copyright 2010 INRA - Cepia Software Platform.


%% Constructors
methods
    function this = CenteredMotionTransform2D(varargin)
        % Create a new model for 2D motion transform model
        % (defined by 1 rotation angle and 2 translation parameters)

        % call parent constructor for initializing center
        this = this@CenteredTransformAbstract([0 0]);
        
        % call parent constructor for initializing parameters
        this = this@ParametricTransform([0 0 0]);
        
        % setup default parameters
        this.params = [0 0 0];
        this.center = [0 0];

        if ~isempty(varargin)
            % extract first argument, and try to interpret
            var = varargin{1};
            if isa(var, 'CenteredMotionTransform2D')
                this.params = var.params;
                
            elseif isnumeric(var)
                if length(var) == 1 && var(1) ~= 2
                    % if argument is a scalar, it corresponds to dimension
                    error('CenteredMotionTransform2D is defined only for 2 dimensions');
                    
                elseif length(var) == 3
                    % if argument is a row vector, it is the parameters
                    this.params = var(1:3);
                else
                    error('Please specify angle in degrees and vector');
                end
            else
                error('Unable to understand input arguments');
            end
        end
        
        % eventually parse additional arguments
        if nargin > 2
            if strcmp(varargin{2}, 'center')
                % setup rotation center
                this.center = varargin{3};
            end
        end
        
        % setup parameter names
        this.paramNames = {'Theta (°)', 'X shift', 'Y shift'};
                
    end % constructor declaration
end

%% Standard methods
methods        
    function dim = getDimension(this) %#ok<MANU>
        dim = 2;
    end

    function mat = affineMatrix(this)
        % Compute affine matrix associated with this transform
        
        % converts to radians
        theta = this.params(1) * pi / 180;

        % pre-computations
        cot = cos(theta);
        sit = sin(theta);

        % center coordinates
        cx = this.center(1);
        cy = this.center(2);

        % translation vector
        ux = this.params(2);
        uy = this.params(3);
        
        % build matrix
        mat = [...
        	cot -sit ( cx*(1-cot) + cy*sit     + ux) ; ...
        	sit +cot (-cx*sit     + cy*(1-cot) + uy) ; ...
            0 0 1 ];

    end
  
    function jacobian = parametricJacobian(this, x, varargin)
        % Compute jacobian matrix, i.e. derivatives for each parameter
       
        % extract coordinate of input point(s)
        if isempty(varargin)
            y = x(:,2);
            x = x(:,1);
        else
            y = varargin{1};
        end
        
        % convert angles to degrees
        theta = this.params(1) * pi / 180;
        
        % precompute angle functions
        cot = cos(theta);
        sit = sin(theta);
        
        % recenter the point coordinates
        x = x - this.center(1);
        y = y - this.center(2);
        
        % compute parametric jacobian
        jacobian = [...
            (-sit*x - cot*y) 1 0 ;
            ( cot*x - sit*y) 0 1 ];
    end

end % methods

%% I/O Methods
methods (Static)
    function transfo = readFromFile(fileName)
        % Read transform from the given file name.
        % Returns a new instance of CenteredMotionTransform2D.
        %
        % Example
        %   TRANSFO = CenteredMotionTransform2D.readFromFile('transfo.txt');
        
        map = readPropertyFile(fileName);
        transfo = CenteredMotionTransform2D.createFromPropertyMap(map);
    end
    
    function transfo = createFromPropertyMap(map)
        % Create a new transform from a set of properties
        
        transfo = CenteredMotionTransform2D();
        
        nbParams = str2double(map('TransformParameterNumber'));
        
        trParams = map('TransformParameters');
        trParams= cellfun(@str2double, regexp(trParams, '\s*', 'split'));
        
        if nbParams ~= length(trParams)
            error('Wrong number of parameters');
        end
        
        setParameters(transfo, trParams);
        
        center = map('TransformCenter');
        center = cellfun(@str2double, regexp(center, '\s*', 'split'));
        transfo.center = center;
    end
end


%% Serialization methods
methods
    function str = toStruct(this)
        % Converts to a structure to facilitate serialization
        str = struct('type', 'CenteredMotionTransform2D', ...
            'center', this.center, ...
            'parameters', this.params);
        
    end
end
methods (Static)
    function transfo = fromStruct(str)
        % Creates a new instance from a structure
        params = str.parameters;
        transfo = CenteredMotionTransform2D(params, 'center', str.center);
    end
end

end % classdef