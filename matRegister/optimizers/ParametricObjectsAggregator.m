classdef ParametricObjectsAggregator < ParametricObject
% Concatenates several parametric objects into a single one.
%
%   output = ParametricObjectsAggregator(input)
%
%   Example
%   ParametricObjectsAggregator
%
%   See also
%

% ------
% Author: David Legland
% e-mail: david.legland@inra.fr
% Created: 2010-11-03,    using Matlab 7.9.0.529 (R2009b)
% Copyright 2010 INRA - Cepia Software Platform.

%% Properties
properties
    % The set of inner parametric objects
    Parametrics;
    
end
 
%% Constructor
methods
    function obj = ParametricObjectsAggregator(varargin)
        % usage:
        % P = ParametricObjectsAggregator(PARAMETRICS);
        
        if nargin ~= 1
            error('Need one input argument');
        end
        
        var = varargin{1};
        if ~iscell(var)
            error('First argument must be a cell array of Parametric objects');
        end
        
        obj.Parametrics = var;
        
    end % constructor
    
end % construction function


%% Provate static methods for inner computations
methods (Static, Access=private)
    function name = formatParameterName(name, i)
        name = sprintf('%s (%d)', name, i);
    end
end


%% Methods implementing the parametric Object interface
methods
    function params = getParameters(obj)
        % Concatenate all parameters vectors into one
        
        nObj = length(obj.Parametrics);
        
        % compute the length of final parameter array
        nOP = getParameterLength(obj);
        
        % concatenate all parameters
        params = zeros(1, nOP);
        ind = 0;
        for i = 1:nObj
            param_i = getParameters(obj.Parametrics{i});
            nP = length(param_i);
            params(ind+(1:nP)) = param_i;
            ind = ind + nP;
        end
        
    end
    
    function setParameters(obj, params)
        % Changes the parameter vector of the transform
        
        nObj = length(obj.Parametrics);
        
        % dispatch parameters to children parametric objects
        ind = 0;
        for i=1:nObj
            item = obj.Parametrics{i};
            nP = getParameterLength(item);
            param_i = params(ind+(1:nP));
            setParameters(item, param_i);            
            ind = ind + nP;
        end
    end
    
    function nOP = getParameterLength(obj)
        % Returns the total length of the parameter array
        
        nObj = length(obj.Parametrics);
        
        % compute the length of final parameter array
        nOP = 0;
        for i = 1:nObj
            nOP = nOP + getParameterLength(obj.Parametrics{i});
        end
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
        
        % iterate over parametric objects
        nOP = 0;
        for i=1:nObj
            item = obj.Parametrics{i};
            nP = getParameterLength(item);
            
            if (nOP+nP) >= paramIndex
                name = getParameterName(item, paramIndex-nOP);
                break;
            end
            
            nOP = nOP + nP;
        end
        name = ParametricObjectsAggregator.formatParameterName(name, i);
    end
    
    function names = getParameterNames(obj)
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

        nObj = length(obj.Parametrics);
        
        % iterate over parametric objects
        nOP = getParameterLength(obj);
        
        % concatenate all parameter names
        names = cell(1, nOP);
        ind = 0;
        for i = 1:nObj
            names_i = getParameterNames(obj.Parametrics{i});
            nP = length(names_i);
            for j = 1:nP
                names{ind+j} = ParametricObjectsAggregator.formatParameterName(names_i{j}, i);
            end
            ind = ind + nP;
        end
     end
    
end % general methods
 

end
