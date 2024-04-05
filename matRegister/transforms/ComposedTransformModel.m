classdef ComposedTransformModel < ParametricTransform
% Compose several transforms, the last one being parametric.
%
%   TC = ComposedTransformModel(T0, T1)
%   T0 is the initial transform, T1 is the transform to optimize. The
%   result TC is the combination of the two transforms.
%
%   Example
%   ParametricTransform2D
%
%   See also
%

% ------
% Author: David Legland
% e-mail: david.legland@inra.fr
% Created: 2010-04-09,    using Matlab 7.9.0.529 (R2009b)
% Copyright 2010 INRA - Cepia Software Platform.

properties
    % the set of encapsulated transforms
    Transforms;
end

%% Constructor
methods
    function obj = ComposedTransformModel(varargin)
        % class constructor.
        % T = ComposedTransformModel(T1, T2);
        
        if isa(varargin{1}, 'ComposedTransform') && nargin == 1
            % copy constructor
            var = varargin{1};
            obj.Transforms = var.Transforms;
            
        elseif isa(varargin{1}, 'Transform')
            % initialize tranform array
            nbTrans = length(varargin);
            obj.Transforms = cell(nbTrans, 1);
            for i = 1:nbTrans
                obj.Transforms{i} = varargin{i};
            end
            
        else
            error('Wrong parameter when constructing a Composed transform');
        end
        
        if ~isa(obj.Transforms{end}, 'ParametricTransform')
            error('Last transform must be a ParametricTransform');
        end
        
    end % constructor
end

%% Implements abstract methods from ParametricTransform
methods
    function jacobian = parametricJacobian(obj, point, varargin)
        % Compute jacobian matrix, i.e. derivatives for coordinate.
        % jacob(i,j) = d x_i / d x_j

        nTransfos = length(obj.Transforms);
        
        % first, transform points
        for i = 1:nTransfos-1
            point = TransformPoint(obj.Transforms{i}, point, varargin{:});
        end

        % then, compute parametric jacobian of the last transform
        jacobian = getParametricJacobian(obj.Transforms{end}, point, varargin{:});
    end

    function transfo = clone(obj)
        transfoList = obj.Transforms;
        transfoList(end) = clone(transfoList(end));
        transfo = ComposedTransformModel(transfoList);
    end
end


%% Overrides several methods from ParametricTransform 

% The goal is to to manipulate params of the last transform instead of
% local parameters

methods
    function p = getParameters(obj)
        % Returns the parameter vector of the transform.
        p = obj.Transforms{end}.Params;
    end
    
    function setParameters(obj, params)
        % Changes the parameter vector of the transform.
        obj.Transforms{end}.Params = params;
    end
    
    function Np = getParameterLength(obj)
        % Returns the length of the vector parameter.
        Np = length(obj.Transforms{end}.Params);
    end
    
    function name = getParameterName(obj, paramIndex)
        % Return the name of the i-th parameter.
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
        if paramIndex > length(obj.Transforms{end}.Params)
            error('Index greater than the number of parameters');
        end
        
        % return a parameter name if it was initialized
        name = '';
        if paramIndex <= length(obj.Transforms{end}.ParamNames)
            name = obj.Transforms{end}.ParamNames{paramIndex};
        end
    end
    
    function name = getParameterNames(obj)
        % Return the names of all parameters in a cell array of strings.
        %
        % NAMES = Transfo.getParameterNames();
        % 
        % Example:
        % T = TranslationModel([10 20]);
        % names = T.getParameterNames()
        % ans = 
        %   'X shift'   'Y shift'
        %
        
        name = obj.Transforms{end}.ParamNames;
    end
end % overridden methods

%% Methods implementing Transform interface
methods
    function dim = getDimension(obj)
        dim = getDimension(obj.Transforms{1});
    end

    function point = transformPoint(obj, point)
        for i = 1:length(obj.Transforms)
            point = TransformPoint(obj.Transforms{i}, point);
        end
    end
    
    function vector = transformVector(obj, vector, position)
        for i = 1:length(obj.Transforms)
            vector = TransformVector(obj.Transforms{i}, vector, position);
        end
    end
    
    function jacobian = jacobianMatrix(obj, point, varargin)
        % Compute jacobian matrix, i.e. derivatives for coordinate.
        % jacob(i,j) = d x_i / d x_j
        
        jacobian = jacobianMatrix(obj.Transforms{1}, point);
        for i = 2:length(obj.Transforms)
            jacobian = jacobianMatrix(obj.Transforms{i}, point, varargin{:}) * jacobian;
        end
    end
end % methods


end % classdef
