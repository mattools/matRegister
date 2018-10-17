classdef ParametricTransform < Transform & ParametricObject
%PARAMETRICTRANSFORM  Abstract class for parametric transform ND->ND
%
%   Superclass for transformation model that can be represented by a vector
%   of parameters, and used as input in optimization procedures.
%

% ------
% Author: David Legland
% e-mail: david.legland@inra.fr
% Created: 2010-04-09,    using Matlab 7.9.0.529 (R2009b)
% Copyright 2010 INRA - Cepia Software Platform.

%% Properties definition
properties
    % the set of inner parameters of the transform
    params;
    
    % the name of each parameter, stored to facilitate automatic plotting.
    % Given as row vector of char arrays.
    paramNames = {};
end

%% Constructor (protected)
methods (Access = protected)
    function this = ParametricTransform(varargin)
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
                this.params = zeros(1, var);
            else
                this.params = var;
            end
        end
        
        % parameter names ciould also be specified
        if nargin > 1
            this.paramNames = varargin{2};
        end
    end
end

%% Methods for managing parameters
methods
    function p = getParameters(this)
        % Returns the parameter vector of the transform
        p = this.params;
    end
    
    function setParameters(this, params)
        % Changes the parameter vector of the transform
        this.params = params;
    end
    
    function Np = getParameterLength(this)
        % Returns the length of the vector parameter
        Np = length(this.params);
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
        if paramIndex > length(this.params)
            error('Index greater than the number of parameters');
        end
        
        % return a parameter name if it was initialized
        name = '';
        if paramIndex <= length(this.paramNames)
            name = this.paramNames{paramIndex};
        end
    end
    
    function name = getParameterNames(this)
        % Return the names of all parameters in a cell array of strings
        %
        % NAMES = Transfo.getParameterNames();
        % 
        % Example:
        % T = TranslationModel([10 20]);
        % names = T.getParameterNames()
        % ans = 
        %   'X shift'   'Y shift'
        %
        
        name = this.paramNames;
    end
    
    function jac = getParametricJacobian(this, x, varargin)
        warning('deprecated: use parametricJacobian method instead');
        jac = parametricJacobian(this, x, varargin{:});
    end
end % methods

%% Abstract methods
methods (Abstract)
    parametricJacobian(this, x, varargin)
    % Compute jacobian matrix, i.e. derivatives for each parameter
    
end % abstract methods 


%% I/O Methods
methods
    function writeToFile(this, file)
        % Write transform parameter to the given file handle
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
        
        nDims = getDimension(this);
        
        fprintf(file, 'TransformType = %s\n', class(this));
        fprintf(file, 'TransformDimension = %d\n', nDims);
        
        nParams = length(this.params);
        fprintf(file, 'TransformParameterNumber = %d \n', nParams);
        
        pattern = ['TransformParameters =', repmat(' %g', 1, nParams) '\n'];
        fprintf(file, pattern, this.params);
        
        % close file
        if closeFile
            fclose(file);
        end
    end
end

end % classdef
