classdef OptimizedValueEvolutionDisplay < OptimizationListener
% Display the evolution of optimized value.
%
%   output = OptimizedValueEvolutionDisplay(input)
%
%   Example
%   OptimizedValueEvolutionDisplay
%
%   See also
%     OptimizationListener, ParametricFunctionEvolutionDisplay
%     ParametersEvolutionDisplay 
%

% ------
% Author: David Legland
% e-mail: david.legland@inra.fr
% Created: 2010-10-26,    using Matlab 7.9.0.529 (R2009b)
% Copyright 2010 INRA - Cepia Software Platform.

properties
    AxisHandle;
    
    ValueArray;
    
    AxisTitle = '';
end

methods
    function obj = OptimizedValueEvolutionDisplay(varargin)
        % Create a new OptimizedValueEvolutionDisplay
        %
        % Example
        %   valueDisplay = OptimizedValueEvolutionDisplay(gcf, 'Image Metric');
        %   addOptimizationListener(optimizer, valueDisplay);

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
        
        if nargin > 1
            obj.AxisTitle = varargin{2};
        end
        
    end % end constructor
    
end

methods
    function optimizationStarted(obj, src, event) %#ok<*INUSD>
        
        % Initialize the parameter array
        obj.ValueArray = src.Value;
    end
    
    function optimizationIterated(obj, src, event)
        
        
        % append current value to the value array
        obj.ValueArray = [obj.ValueArray; src.Value];
        
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
