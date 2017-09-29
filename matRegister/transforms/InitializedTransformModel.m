classdef InitializedTransformModel < ParametricTransform
% Encapsulation of a parametric and an initial transform
%
%   TC = InitializedTransformModel(T0, T1)
%   T0 is the initial transform, T1 is the transform to optimize. The
%   result TC is the combination of the two transforms.
%
%   Example
%   InitializedTransformModel
%
%   See also
%

% ------
% Author: David Legland
% e-mail: david.legland@inra.fr
% Created: 2010-04-09,    using Matlab 7.9.0.529 (R2009b)
% Copyright 2010 INRA - Cepia Software Platform.

properties
    % the initial transform
    initial;
    
    % the parametric transform to optimize
    transform;
end

%% Constructor
methods
    function this = InitializedTransformModel(varargin)
        % class constructor
        % T = InitializedTransformModel(T1, T2);
        
        if isa(varargin{1}, 'ComposedTransform') && nargin == 1
            % copy constructor
            var = varargin{1};
            this.transforms = var.transforms;
            
        elseif nargin == 2
            % initialization constructor
            this.initial = varargin{1};
            this.transform = varargin{2};
            
        else
            error('Wrong parameter when constructing a Composed transform');
        end
        
        if ~isa(this.transform, 'ParametricTransform')
            error('Last transform must be a ParametricTransform');
        end
        
    end % constructor
end


%% Methods implementing Transform interface
methods
    function dim = getDimension(this)
        dim = getDimension(this.initial);
    end

    function point = transformPoint(this, point)
        point = transformPoint(this.transform, transformPoint(this.initial, point));
    end
    
    function vector = transformVector(this, vector, position)
        error('Unimplemened method');
    end
    
    function jacobian = getJacobian(this, point, varargin)
        % Compute jacobian matrix, i.e. derivatives for coordinate
        % jacob(i,j) = d x_i / d x_j
        
        jacobian = this.initial.getJacobian(point);
        jacobian = this.transform.getJacobian(point) * jacobian;
    end
    
    function jacobian = getParametricJacobian(this, point, varargin)
        % Compute jacobian matrix, i.e. derivatives for coordinate
        % jacob(i,j) = d x_i / d x_j

        % first, transform points
        point = transformPoint(this.initial, point, varargin{:});
        
        % then, compute parametric jacobian of the last transform
        jacobian = getParametricJacobian(this.transform, point, varargin{:});
        
    end
end % methods


%% Overrides several methods from ParametricTransform 

% The goal is to to manipulate params of the last transform instead of
% local parameters

methods
    function p = getParameters(this)
        % Returns the parameter vector of the transform
        p = this.transform.params;
    end
    
    function setParameters(this, params)
        % Changes the parameter vector of the transform
        this.transform.params = params;
    end
    
    function Np = getParameterLength(this)
        % Returns the length of the vector parameter
        Np = length(this.transform.params);
    end
    
    function name = getParameterName(this, paramIndex)
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
        if paramIndex > length(this.transform.params)
            error('Index greater than the number of parameters');
        end
        
        % return a parameter name if it was initialized
        name = '';
        if paramIndex <= length(this.transform.paramNames)
            name = this.transform.paramNames{paramIndex};
        end
    end
    
    function names = getParameterNames(this)
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
        
        names = this.transform.paramNames;
    end
end % overridden methods

end % classdef
