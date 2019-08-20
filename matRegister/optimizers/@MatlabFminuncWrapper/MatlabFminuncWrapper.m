classdef MatlabFminuncWrapper < Optimizer
%MATLABFMINUNCWRAPPER  One-line description here, please.
%
%   output = MatlabFminuncWrapper(input)
%
%   Example
%   MatlabFminuncWrapper
%
%   See also
%

% ------
% Author: David Legland
% e-mail: david.legland@inra.fr
% Created: 2011-01-12,    using Matlab 7.9.0.529 (R2009b)
% Copyright 2011 INRA - Cepia Software Platform.


%% Properties
properties
   NIters = 200;
end 

%% Constructor
methods
    function obj = MatlabFminuncWrapper(varargin)
        
        obj = obj@Optimizer();
    end % constructor 

end % construction function

end % classdef
