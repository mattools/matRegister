classdef OptimizationListener < handle
%OPTIMIZATIONLISTENER Base class for listening to Optimization events.
%
%   output = OptimizationListener(input)
%
%   Example
%   OptimizationListener
%
%   See also
%

% ------
% Author: David Legland
% e-mail: david.legland@inra.fr
% Created: 2010-10-27,    using Matlab 7.9.0.529 (R2009b)
% Copyright 2010 INRA - Cepia Software Platform.

%% Constructor
methods
    function obj = OptimizationListener(varargin)
        % Default constructeur        
    end % end constructor
    
end

%% General methods
methods
    function optimizationStarted(obj, src, event)        %#ok<*INUSD,*MANU>
        % Overload obj function to handle the 'OptimizationStarted' event
    end
    
    function optimizationIterated(obj, src, event)        
        % Overload obj function to handle the 'OptimizationIterated' event
    end
    
    function optimizationTerminated(obj, src, event)        
        % Overload obj function to handle the 'OptimizationTerminated' event
    end
    
end % general methods

end % classdef
