classdef CenteredEulerTransform3D < AffineTransform & ParametricTransform & CenteredTransformAbstract
%Transformation model for a centered 3D rotation followed by a translation.
%   
%   Inner optimisable parameters of the transform have the following form:
%   params[1] = phi, rotation angle around Ox, in degrees
%   params[2] = theta, rotation angle around Oy, in degrees
%   params[3] = psi, rotation angle around Oz, in degrees
%   params[4] = tx
%   params[5] = ty
%   params[6] = tz
% 
%   Creation:
%   Trans = CenteredEulerTransform3D();
%   

% ------
% Author: David Legland
% e-mail: david.legland@inra.fr
% Created: 2010-07-27,    using Matlab 7.9.0.529 (R2009b)
% Copyright 2010 INRA - Cepia Software Platform.


%% Constructors
methods
    function obj = CenteredEulerTransform3D(varargin)
        % Create a new centered motion transform
        
        obj = obj@CenteredTransformAbstract([0 0 0]);
        
        obj.Params = zeros(1, 6);
        
        if ~isempty(varargin)
            % extract first argument, and try to interpret
            var = varargin{1};
            if isa(var, 'CenteredEulerTransform3D')
                % copy constructor
                obj.Params = var.Params;
                
            elseif isnumeric(var)
                len = length(var);
                if len == 6
                    % all parameters are specified
                    obj.Params = var;
                elseif len == 3
                    % specify only angles, initialize translation to 0
                    obj.Params = [var 0 0 0];
                elseif len == 1
                    % give the possibility to call constructor with the
                    % dimension of image
                    if var ~= 3
                        error('Defined only for 3 dimensions');
                    end
                else
                    error('Please specify 3 angles and a 3D vector');
                end
                
            else
                error('Unable to understand input arguments');
            end
        end
        
        % eventually parse additional arguments
        if nargin > 2
            if strcmp(varargin{2}, 'Center')
                % setup rotation center
                obj.Center = varargin{3};
            end
        end
        
        % setup parameter names
        obj.ParamNames = {...
            'Rot X (°)', 'Rot Y (°)', 'Rot Z (°)', ...
            'X shift', 'Y shift', 'Z shift'};
                
    end % constructor declaration
end

%% Standard methods
methods    
    function dim = getDimension(obj) %#ok<MANU>
        dim = 3;
    end

    function initFromTranslation(obj, vector)
        % Initialize parameters from a translation vector
        %
        % Example
        % T = CenteredEulerTransform3D();
        % T.initFromTranslation([10 15 20]);
        % T.getParameters()
        % ans = 
        %     [0 0 0 10 15 20]
        %
        obj.Parameters = [0 0 0 vector];
    end
    
    function mat = affineMatrix(obj)
        % Compute affine matrix associated with obj transform
        
        % extract angles and convert to radians
        deg2rad = pi / 180;
        phi     = obj.Params(1) * deg2rad;
        theta   = obj.Params(2) * deg2rad;
        psi     = obj.Params(3) * deg2rad;

%         % compute compound rotation matrix around center
%         Rx = createRotationOx(phi);
%         Ry = createRotationOy(theta);
%         Rz = createRotationOz(psi);
%         mat = recenterTransform3d(Rz * Ry * Rx, obj.Center);
%         % add translation
%         mat = createTranslation3d(obj.Params(4:6)) * mat;

        % pre-computations of trigonometric functions
        a = cos(phi);      b = sin(phi);
        c = cos(theta);    d = sin(theta);
        e = cos(psi);      f = sin(psi);
        
        % base matrix
        bd = b * d;
        ad = a * d;
        mat = [c*e bd*e-a*f ad*e+b*f; c*f bd*f+a*e ad*f-b*e; -d b*c a*c];
        
        p0 = obj.Center';
        shift = p0 - mat*p0 + obj.Params(4:6)';
        mat = [mat shift ; 0 0 0 1];
    end
  
    function jacobian = parametricJacobian(obj, x, varargin)
        % Compute jacobian matrix, i.e. derivatives for each parameter
       
        % extract coordinate of input point(s)
        if isempty(varargin)
            y = x(:,2);
            z = x(:,3);
            x = x(:,1);
        else
            y = varargin{1};
            z = varargin{2};
        end
        
        % extract angles in degrees
        k = pi / 180;
        phi     = obj.Params(1) * k;
        theta   = obj.Params(2) * k;
        psi     = obj.Params(3) * k;

        % pre-computations of trigonometric functions (in degrees)
        cx = cos(phi);      sx = sin(phi);
        cy = cos(theta);    sy = sin(theta);
        cz = cos(psi);      sz = sin(psi);

        % jacobians are computed with respect to transformation center
        x = x - obj.Center(1);
        y = y - obj.Center(2);
        z = z - obj.Center(3);
        
        % compute the Jacobian matrix using pre-computed elements
        jacobian = [[...
            ((cz*sy*cx + sz*sx)*y + (sz*cx - cz*sy*sx)*z), ...
            ((sz*sy*cx - cz*sx)*y - (sz*sy*sx + cz*cx)*z), ...
            (             cy*cx*y -              cy*sx*z); ...
            ...
            (-cz*sy*x + cz*cy*sx*y + cz*cy*cx*z), ...
            (-sz*sy*x + sz*cy*sx*y + sz*cy*cx*z), ...
            (-cy*x    -    sy*sx*y -    sy*cx*z); ...
            ...
            (-sz*cy*x - (sz*sy*sx + cz*cx)*y + (cz*sx - sz*sy*cx)*z), ...
            ( cz*cy*x + (cz*sy*sx - sz*cx)*y + (cz*sy*cx + sz*sx)*z), ...
            0; ...
            ]' ...
            eye(3, 3)]; % correspond to translation part of the transform
    end

end % methods

%% I/O Methods
methods (Static)
    function transfo = readFromFile(fileName)
        % Read transform from the given file name.
        % Returns a new instance of CenteredEulerTransform3D.
        %
        % Example
        %   TRANSFO = CenteredEulerTransform3D.readFromFile('transfo.txt');
        
        map = readPropertyFile(fileName);
        transfo = CenteredEulerTransform3D.createFromPropertyMap(map);
    end
    
    function transfo = createFromPropertyMap(map)
        % Create a new transform from a set of properties
        
        transfo = CenteredEulerTransform3D();
        
        nbParams = str2double(map('TransformParameterNumber'));
        
        trParams = map('TransformParameters');
        trParams= cellfun(@str2double, regexp(trParams, '\s*', 'split'));
        
        if nbParams ~= length(trParams)
            error('Wrong number of parameters');
        end
        
        setParameters(transfo, trParams);
        
        center = map('TransformCenter');
        center = cellfun(@str2double, regexp(center, '\s*', 'split'));
        transfo.Center = center;
    end
end


%% Serialization methods
methods
    function str = toStruct(obj)
        % Converts to a structure to facilitate serialization
        str = struct('Type', 'CenteredEulerTransform3D', ...
            'Center', obj.Center, ...
            'Parameters', obj.Params);
        
    end
end
methods (Static)
    function transfo = fromStruct(str)
        % Creates a new instance from a structure
        params = str.Parameters;
        transfo = CenteredEulerTransform3D(params, 'Center', str.Center);
    end
end

end % classdef