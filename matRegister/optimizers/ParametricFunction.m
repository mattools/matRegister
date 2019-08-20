classdef ParametricFunction < handle
% Abstract classe for parametric objects that compute a value from params.
%
%   output = ParametricFunction(input)
%
%   Example
%   ParametricFunction
%
%   See also
%

% ------
% Author: David Legland
% e-mail: david.legland@inra.fr
% Created: 2010-10-27,    using Matlab 7.9.0.529 (R2009b)
% Copyright 2010 INRA - Cepia Software Platform.


%% Abstract methods
methods (Abstract)
    setParameters(obj, params)
    % Setup the parameters of the object
  
%     params = getParameters(obj)
%     % Returns the parameter vector associated to obj transform
    
    res = computeValue(obj)
    % Compute the value using the actual parameters (do not change state)
    
end % abstract methods

methods
    function [res, grad] = evaluate(obj, params)
        % basic implementation of evaluate function
        % obj make possible the call in an Optimization procedure.
        setParameters(obj, params);
        if nargout <= 1
            res = computeValue(obj);
        else
            [res, grad] = computeValueAndGradient(obj);
        end
    end
    
end % base implementation

end % classdef


