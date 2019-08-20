classdef ParametersEvolutionDisplay < OptimizationListener
% Display the evolution of optimization parameters.
%
%   output = ParametersEvolutionDisplay(input)
%
%   Example
%   ParametersEvolutionDisplay
%
%   See also
%     OptimizationListener, ParametricFunctionEvolutionDisplay
%     ParametersEvolutionDisplay

% ------
% Author: David Legland
% e-mail: david.legland@inra.fr
% Created: 2010-10-26,    using Matlab 7.9.0.529 (R2009b)
% Copyright 2010 INRA - Cepia Software Platform.

properties
    FigureHandle;
    
    PlotMatrix;
    
    ParamValues;
    
    Labels;
end

methods
    function obj = ParametersEvolutionDisplay(varargin)
        if nargin == 0
            obj.FigureHandle = gcf;
        else
            obj.FigureHandle = varargin{1};
        end
        
        if nargin > 1
            obj.PlotMatrix = varargin{2};
        else
            obj.PlotMatrix = [2 2];
        end
        
        if nargin > 2
            var = varargin{3};
            if isnumeric(var)
                obj.Labels = strtrim(cellstr(num2str((1:10)', 'Param %d')));
            elseif iscell(var)
                obj.Labels = var;
            else
                error('Can not process labels parameter');
            end
            
        else
            obj.Labels = '';
        end
        
    end % end constructor
    
end

methods
    function optimizationStarted(obj, src, event) %#ok<*INUSD>
        % Initialize the parameter array
        params = getParameters(src);
        obj.ParamValues = params;
    end
    
    function optimizationIterated(obj, src, event)
        
        % append current parameters to the parameter array
        params = src.getParameters();
        obj.ParamValues = [obj.ParamValues ; params];
        
        figure(obj.FigureHandle);
        nRows = obj.PlotMatrix(1);
        nCols = obj.PlotMatrix(2);
        nv = size(obj.ParamValues, 1);
        
        for i = 1:length(params)
            subplot(nRows, nCols, i);
            plot(1:nv, obj.ParamValues(:, i));
            xlim([0 nv]);
            
            if ~isempty(obj.Labels)
                title(obj.Labels{i});
            end
        end
        drawnow;
    end
    
    function optimizationTerminated(obj, src, event)
    end
    
end

end
