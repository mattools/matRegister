classdef BaseFunction < handle
% Parent class for objects that can compute a scalar value
%
%   output = BaseFunction(input)
%
%   Example
%   BaseFunction
%
%   See also
%

% ------
% Author: David Legland
% e-mail: david.legland@inra.fr
% Created: 2010-11-03,    using Matlab 7.9.0.529 (R2009b)
% Copyright 2010 INRA - Cepia Software Platform.

methods (Abstract)
    value = computeValue(varargin)
end

end
