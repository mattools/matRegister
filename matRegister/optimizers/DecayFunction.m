classdef DecayFunction < handle
% Abstract class that models decay function for optimizers.
%
%   Class DecayFunction
%
%   Example
%   Decayfunction
%
%   See also
%     ExponentialDecayFunction, ElastixDecayFunction
%

% ------
% Author: David Legland
% e-mail: david.legland@inra.fr
% Created: 2011-08-04,    using Matlab 7.9.0.529 (R2009b)
% Copyright 2011 INRA - Cepia Software Platform.


%% Constructor
methods
    function obj = DecayFunction(varargin)
    % Constructor for Decayfunction class

    end

end % end constructors


%% Methods
methods (Abstract)
    alpha = evaluate(iter)
end % end methods

end % end classdef

