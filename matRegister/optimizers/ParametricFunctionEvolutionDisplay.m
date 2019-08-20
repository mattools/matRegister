classdef ParametricFunctionEvolutionDisplay < OptimizationListener
%PARAMETRICFUNCTIONEVOLUTIONDISPLAY Displays evolution of a parametric function
%
%   output = ParametricFunctionEvolutionDisplay(input)
%
%   Example
%   ParametricFunctionEvolutionDisplay
%
%   See also
%

% ------
% Author: David Legland
% e-mail: david.legland@inra.fr
% Created: 2010-10-26,    using Matlab 7.9.0.529 (R2009b)
% Copyright 2010 INRA - Cepia Software Platform.

properties
    AxisHandle;
    
    ParametricFunction;
    
    ValueArray;
    
    AxisTitle = '';
end

methods
    function obj = ParametricFunctionEvolutionDisplay(varargin)
        
        if nargin < 2
            error('Need at least two input arguments');
        end
        
        % Initialize axis handle
        if ~isempty(varargin)
            var = varargin{1};
            if ~ishandle(var)
                error('First argument must be an axes handle');
            end
            
            if strcmp(get(var, 'Type'), 'axes')
                obj.AxisHandle = var;
            else
                obj.AxisHandle = gca;
            end
        else
            obj.AxisHandle = gcf;
        end

        obj.ParametricFunction = varargin{2};
        
        if nargin > 2
            obj.AxisTitle = varargin{3};
        end
        
    end % end constructor
    
end

methods
    function optimizationStarted(obj, src, event) %#ok<*INUSD>
        
        % Initialize the parameter array
        obj.ValueArray = src.Value;
    end
    
    function optimizationIterated(obj, src, event)
        
        
        % compute current value of the function
        value = computeValue(obj.ParametricFunction);
        
        % append current value to the value array
        obj.ValueArray = [obj.ValueArray; value];
        
        % display current list of values
        nv = length(obj.ValueArray);
        plot(obj.AxisHandle, 1:nv, obj.ValueArray);
        set(obj.AxisHandle, 'xlim', [0 nv]);
        
        % decorate
        if ~isempty(obj.AxisTitle)
            title(obj.AxisHandle, obj.AxisTitle);
        end
        
        % refresh display
        drawnow expose;
    end
    
end

end
