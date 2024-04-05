classdef MotionModel2D < AffineTransform & ParametricTransform
%Transformation model for a centered rotation followed by a translation.
%   
%   Inner optimisable parameters of the transform have the following form:
%   params[1] = tx
%   params[2] = ty
%   params[3] = theta, in degrees
% 
%   TRANS = MotionModel2D();
%   Create using default parameters (all zero).
%
%   TRANS = MotionModel2D(PARAMS);
%   Create a new transform model by initializing parameters. PARAMS is a
%   1-by-3 row vector containing TX, TY, and THETA parameters.
%
%

% ------
% Author: David Legland
% e-mail: david.legland@inra.fr
% Created: 2017-08-04,    using Matlab 7.9.0.529 (R2009b)
% Copyright 2010 INRA - Cepia Software Platform.


%% Constructors
methods
    function obj = MotionModel2D(varargin)
        % Create a new model for 2D motion transform model.
        % (defined by 2 translation parameters and 1 rotation angle)

        % call parent constructor for initializing parameters
        obj = obj@ParametricTransform([0 0 0]);
        
        % setup default parameters
        obj.Params = [0 0 0];

        if ~isempty(varargin)
            % extract first argument, and try to interpret
            var = varargin{1};
            if isa(var, 'MotionModel2D')
                % copy constructor
                obj.Params = var.Params;
                
            elseif isnumeric(var)
                if isscalar(var) && var(1) ~= 2
                    % if argument is a scalar, it corresponds to dimension
                    error('MotionModel2D is defined only for 2 dimensions');
                    
                elseif length(var) == 3
                    % if argument is a row vector, it is the parameters
                    obj.Params = var(1:3);
                else
                    error('Please specify angle in degrees and vector');
                end
            else
                error('Unable to understand input arguments');
            end
        end
       
        % setup parameter names
        obj.ParamNames = {'X shift', 'Y shift', 'Theta (�)'};
                
    end % constructor declaration
end


%% Implementation of ParametricTransform methods
methods
    function jacobian = parametricJacobian(obj, x, varargin)
        % Compute jacobian matrix, i.e. derivatives for each parameter.
       
        % extract coordinate of input point(s)
        if isempty(varargin)
            y = x(:,2);
            x = x(:,1);
        else
            y = varargin{1};
        end
        
        % convert angles to radians
        theta = obj.Params(3) * pi / 180;
        
        % precompute angle functions
        cot = cos(theta);
        sit = sin(theta);
        
        % compute parametric jacobian
        jacobian = [...
            1 0 (-sit*x - cot*y) ;
            0 1 ( cot*x - sit*y) ];
    end

    function transfo = clone(obj)
        transfo = MotionModel2D(obj.Params);
    end
end


%% Standard methods
methods        
    function mat = affineMatrix(obj)
        % Compute affine matrix associated with obj transform.
        
        % translation vector
        tx = obj.Params(1);
        ty = obj.Params(2);
        
        % rotation angle, converted to radians
        theta = obj.Params(3) * pi / 180;

        % pre-computations
        cot = cos(theta);
        sit = sin(theta);

        % build matrix
        mat = [...
        	cot -sit tx ; ...
        	sit +cot ty ; ...
             0    0   1 ];
    end
  
    function dim = getDimension(obj) %#ok<MANU>
        dim = 2;
    end

end % methods

% %% I/O Methods
% methods (Static)
%     function transfo = readFromFile(fileName)
%         % Read transform from the given file name.
%         % Returns a new instance of MotionModel2D.
%         %
%         % Example
%         %   TRANSFO = MotionModel2D.readFromFile('transfo.txt');
%         
%         map = readPropertyFile(fileName);
%         transfo = MotionModel2D.createFromPropertyMap(map);
%     end
%     
%     function transfo = createFromPropertyMap(map)
%         % Create a new transform from a set of properties
%         
%         transfo = MotionModel2D();
%         
%         nbParams = str2double(map('TransformParameterNumber'));
%         
%         trParams = map('TransformParameters');
%         trParams= cellfun(@str2double, regexp(trParams, '\s*', 'split'));
%         
%         if nbParams ~= length(trParams)
%             error('Wrong number of parameters');
%         end
%         
%         setParameters(transfo, trParams);
%         
%         center = map('TransformCenter');
%         center = cellfun(@str2double, regexp(center, '\s*', 'split'));
%         transfo.center = center;
%     end
% end

%% Serialization methods
methods
    function str = toStruct(obj)
        % Converts to a structure to facilitate serialization.
        str = struct('Type', 'MotionModel2D', ...
            'Parameters', obj.Params);
    end
end
methods (Static)
    function motion = fromStruct(str)
        % Creates a new instance from a structure.
        if isfield(str, 'Parameters')
            params = str.Parameters;
        elseif isfield(str, 'parameters')
            params = str.parameters;
        end
        motion = MotionModel2D(params);
    end
end

end % classdef