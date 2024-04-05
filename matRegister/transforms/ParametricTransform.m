classdef ParametricTransform < Transform & ParametricObject
%PARAMETRICTRANSFORM  Abstract class for parametric transform ND->ND.
%
%   Superclass for transformation model that can be represented by a vector
%   of parameters, and used as input in optimization procedures.
%

% ------
% Author: David Legland
% e-mail: david.legland@inrae.fr
% Created: 2010-04-09,    using Matlab 7.9.0.529 (R2009b)
% Copyright 2010 INRA - Cepia Software Platform.


%% Properties definition
properties
    % the set of inner parameters of the transform
    Params;
    
    % the name of each parameter, stored to facilitate automatic plotting.
    % Given as row vector of char arrays.
    ParamNames = {};
end


%% Constructor (protected)
methods (Access = protected)
    function obj = ParametricTransform(varargin)
        % Construct a new parametric transform.
        %
        % TRANSFO = ParametricTransform(NP);
        % Specifies the number of parameters.
        %
        % TRANSFO = ParametricTransform(PARAMS);
        % Specifies the initial parameter vector. PARAMS should by a row
        % vector.
        %
        % TRANSFO = ParametricTransform(..., NAMES);
        % Also specifies the parameter names.
        %
        
        % initialize parameter vector
        if ~isempty(varargin)
            var = varargin{1};
            if isscalar(var)
                obj.Params = zeros(1, var);
            else
                obj.Params = var;
            end
        end
        
        % parameter names ciould also be specified
        if nargin > 1
            obj.ParamNames = varargin{2};
        end
    end
end


%% Methods for managing parameters
methods
    function p = getParameters(obj)
        % Returns the parameter vector of the transform.
        p = obj.Params;
    end
    
    function setParameters(obj, params)
        % Changes the parameter vector of the transform.
        obj.Params = params;
    end
    
    function Np = getParameterLength(obj)
        % Returns the length of the vector parameter.
        Np = length(obj.Params);
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
        if paramIndex > length(obj.Params)
            error('Index greater than the number of parameters');
        end
        
        % return a parameter name if it was initialized
        name = '';
        if paramIndex <= length(obj.ParamNames)
            name = obj.ParamNames{paramIndex};
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
        
        name = obj.ParamNames;
    end
    
    function jac = getParametricJacobian(obj, x, varargin)
        % deprecated: use parametricJacobian instead
        warning('deprecated: use parametricJacobian method instead');
        jac = parametricJacobian(obj, x, varargin{:});
    end
end % methods


%% Abstract methods
methods (Abstract)
    % Compute jacobian matrix, i.e. derivatives for each parameter.
    % 
    % jac = parametricJacobian(transfo, pos);
    % The Jacobian matrix has as many rows as the number of dimensions of
    % the transform (usually 2 or 3), and as many columns as the number of
    % parameters.
    % 
    % Example
    %     % compute parametric jacobian of 2D motion transform model
    %     transfo = MotionModel2D([10 20 30]);
    %     parametricJacobian(transfo, [10 10])    
    %     ans =
    %         1.0000         0  -13.6603
    %              0    1.0000    3.6603
    parametricJacobian(obj, x, varargin)

    % Create a copy of this parametric transform.
    clone(obj);

end % abstract methods 


%% I/O Methods
methods
    function writeToFile(obj, file)
        % Write transform parameter to the given file handle.
        % Assumes file handle is an instance of FileWriter.
        %
        % Example
        %   F = fopen('transfo.txt', 'wt');
        %   fprintf(F, '#--- Transform Parameters ---');
        %   writeToFile(TRANSFO, F);
        %   fclose(F);
        %
        
        closeFile = false;
        if ischar(file)
            file = fopen(file, 'wt');
            closeFile = true;
        end
        
        nDims = getDimension(obj);
        
        fprintf(file, 'TransformType = %s\n', class(obj));
        fprintf(file, 'TransformDimension = %d\n', nDims);
        
        nParams = length(obj.Params);
        fprintf(file, 'TransformParameterNumber = %d \n', nParams);
        
        pattern = ['TransformParameters =', repmat(' %g', 1, nParams) '\n'];
        fprintf(file, pattern, obj.Params);
        
        % close file
        if closeFile
            fclose(file);
        end
    end
end

end % classdef
