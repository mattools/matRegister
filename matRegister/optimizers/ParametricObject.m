classdef ParametricObject < handle
%PARAMETRICOBJECT Base class for all parametric objects
%
%   output = ParametricObject(input)
%
%   Example
%   ParametricObject
%
%   See also
%

% ------
% Author: David Legland
% e-mail: david.legland@inra.fr
% Created: 2010-11-03,    using Matlab 7.9.0.529 (R2009b)
% Copyright 2010 INRA - Cepia Software Platform.


%% Methods for managing parameters
methods (Abstract)
    params = getParameters(obj)
    % Returns the parameter vector of the transform
    
    setParameters(obj, params)
    % Changes the parameter vector of the transform
    
    Np = getParameterLength(obj)
    % Returns the length of the vector parameter
    
    name = getParameterName(obj, paramIndex)
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
    
    name = getParameterNames(obj)
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
    
end % methods

end % classdef
