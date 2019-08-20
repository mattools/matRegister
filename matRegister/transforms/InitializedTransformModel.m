classdef InitializedTransformModel < ParametricTransform
% Encapsulation of a parametric and an initial transform.
%
%   TC = InitializedTransformModel(T0, T1)
%   T0 is the initial transform, T1 is the transform to optimize. The
%   result TC is the combination of the two transforms.
%
%   Example
%   InitializedTransformModel
%
%   See also
%     ParametricTransform, ComposedTransform

% ------
% Author: David Legland
% e-mail: david.legland@inra.fr
% Created: 2010-04-09,    using Matlab 7.9.0.529 (R2009b)
% Copyright 2010 INRA - Cepia Software Platform.

properties
    % the initial transform
    Initial;
    
    % the parametric transform to optimize
    Transform;
end

%% Constructor
methods
    function obj = InitializedTransformModel(varargin)
        % class constructor
        % T = InitializedTransformModel(T1, T2);
        
        if isa(varargin{1}, 'ComposedTransform') && nargin == 1
            % copy constructor
            var = varargin{1};
            obj.Initial = var.Initial;
            obj.Transform = var.Transform;
            
        elseif nargin == 2
            % initialization constructor
            obj.Initial = varargin{1};
            obj.Transform = varargin{2};
            
        else
            error('Wrong parameter when constructing an initialized transform model');
        end
        
        % check validity of arguments
        if ~isa(obj.Initial, 'Transform')
            error('initial transform must be an instance of Transform');
        end
        if ~isa(obj.Transform, 'ParametricTransform')
            error('Last transform must be a ParametricTransform');
        end
        
    end % constructor
end


%% Overrides several methods from ParametricTransform 

% The goal is to to manipulate params of the last transform instead of
% local parameters

methods
    function p = getParameters(obj)
        % Returns the parameter vector of the transform
        p = obj.Transform.Params;
    end
    
    function setParameters(obj, params)
        % Changes the parameter vector of the transform
        obj.Transform.Params = params;
    end
    
    function Np = getParameterLength(obj)
        % Returns the length of the vector parameter
        Np = length(obj.Transform.Params);
    end
    
    function name = getParameterName(obj, paramIndex)
        % Return the name of the i-th parameter
        %
        % NAME = Transfo.getParameterName(PARAM_INDEX);
        % PARAM_INDEX is the parameter index, between 0 and the number of
        % parameters.
        %
        % T = TranslationModel([10 20]);
        % name = T.getParameterName(2);
        % name = 
        %   Y shift
        %
        
        % check index is not too high
        if paramIndex > length(obj.Transform.Params)
            error('Index greater than the number of parameters');
        end
        
        % return a parameter name if it was initialized
        name = '';
        if paramIndex <= length(obj.Transform.ParamNames)
            name = obj.Transform.ParamNames{paramIndex};
        end
    end
    
    function names = getParameterNames(obj)
        % Return the names of all parameters in a cell array of strings
        %
        % NAMES = transfo.getParameterNames();
        % 
        % Example:
        % T = TranslationModel([10 20]);
        % names = T.getParameterNames()
        % ans = 
        %   'X shift'   'Y shift'
        %
        
        names = obj.Transform.ParamNames;
    end
    
    function jacobian = parametricJacobian(obj, point, varargin)
        % Compute jacobian matrix, i.e. derivatives for coordinate
        % jacob(i,j) = d x_i / d x_j

        % first, transform points
        point = transformPoint(obj.Initial, point, varargin{:});
        
        % then, compute parametric jacobian of the last transform
        jacobian = parametricJacobian(obj.Transform, point, varargin{:});
        
    end
end % methods implementing ParametricTransform interface

%% Methods implementing Transform interface
methods
    function point = transformPoint(obj, point)
        point = transformPoint(obj.Transform, transformPoint(obj.Initial, point));
    end
    
    function jacobian = jacobianMatrix(obj, point, varargin)
        % Compute jacobian matrix, i.e. derivatives for coordinate
        % jacob(i,j) = d x_i / d x_j
        
        jacobian = jacobianMatrix(obj.Initial, point);
        jacobian = jacobianMatrix(obj.Transform, point) * jacobian;
    end
    
    function dim = getDimension(obj)
        dim = getDimension(obj.Initial);
    end

end % methods implementing Transform interface


%% Serialization methods

methods
    function str = toStruct(obj)
        % Converts to a structure to facilitate serialization
        str = struct('Type', 'InitializedTransformModel', ...
            'Initial', toStruct(obj.Initial), ...
            'Transform', toStruct(obj.Transform));
    end
end

methods (Static)
    function transfo = fromStruct(str)
        % Creates a new instance from a structure
        init = Transform.fromStruct(str.Initial);
        transform = Transform.fromStruct(str.Transform);
        transfo = InitializedTransformModel(init, transform);
    end
end

end % classdef
