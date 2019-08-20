classdef Optimizer < handle
% General interface for single-valued function optimizers.
%
%   This abstract class is used to define the general contract of
%   optimization algorithms. All implementations of Optimizer abstraction
%   should implement the "startOptimization" method, with two syntaxes:
%   * PARAMS = startOptimization(OPTIM)
%   * [PARAMS, VALUE] = startOptimization(OPTIM)
%
%   It would be nice to also support following syntax:
%   * PARAMS = startOptimization((OPTIM, PARAMS0)
%
%   The Optimizer class implements a list of listeners. 
%   Example of use:
%     fun = @rosenbrock;
%     t0 = [0 0]; deltas = [.01 .01];
%     optimizer = NelderMeadSimplexOptimizer(fun, [0 0], [.01 .01]);
%     listener = OptimizedValueEvolutionDisplay(gca);
%     addOptimizationListener(optimizer, listener);
%     [xOpt value] = optimizer.startOptimization();
%     
%
%   See also
%   optimizers, NelderMeadSimplexOptimizer, MultiLinearSearchOptimizer, 
%   GaussianLinearSearchOptimizer, GradientDescentOptimizer, 
%   BoundedMultiLinearOptimizer, MatlabSimplexOptimizer, MatlabFminuncWrapper
%

% ------
% Author: David Legland
% e-mail: david.legland@inra.fr
% Created: 2010-10-06,    using Matlab 7.9.0.529 (R2009b)
% Copyright 2010 INRA - Cepia Software Platform.

properties
    % the function to minimize
    CostFunction;
    
    % the initial set of parameters
    InitialParameters = [];
    
    % the current set of parameters
    Params;

    % the current value
    Value; 
    
    % Some scaling of the parameters for homogeneization
    % (parameters will be divided by corresponding scale)
    ParameterScales = [];
    
    % the function that will be called at each iteration
    % (the use of Optimization Events is more convenient)
    OutputFunction = [];
    
    % Specifies which information will be displayed at each iteration
    % Can be one of:
    % 'iter'    display information at each iteration
    % 'final'   display only message of convergence
    % 'notify'  display message if convergence failed
    % 'off'     does not display anything
    DisplayMode = 'notify';
end

%% Events
events
    % notified when the optimization starts
    OptimizationStarted
    
    % notified when the an iteration is run
    OptimizationIterated
    
    % notified when the optimization finishes
    OptimizationTerminated
end

%% Constructor
methods (Access = protected)
    function obj = Optimizer(costFun, params0, varargin)
        % Initialize a new Optimizer
        %
        % Usage (in sub-class constructor):
        % obj = obj@Optimizer();
        % obj = obj@Optimizer(FUN, PARAMS);
        % FUN is either a function handle, or an instance of CostFunction
        % PARAMS is the initial set of parameters
        %
        % See Also
        %   setCostFunction, setInitialParameters
        
        if nargin == 0
            return;
        end
        
        setCostFunction(obj, costFun);
        setInitialParameters(obj, params0);
        setParameters(obj, params0);
        
    end 
        
end % Constructors

%% Abstract methods
methods (Abstract)
    varargout = startOptimization(varargin)
    %STARTOPTIMIZATION Start the optimizer and iterate until an end condition is reached
end

%% General methods
methods
    function params = getInitialParameters(obj)
        params = obj.InitialParameters;
    end
    
    function setInitialParameters(obj, params0)
        obj.InitialParameters = params0;
    end
    
    function params = getParameters(obj)
        params = obj.Params;
    end
    
    function setParameters(obj, params)
        obj.Params = params;
    end
    
    function scales = getParameterScales(obj)
        scales = obj.ParameterScales;
    end
    
    function setParameterScales(obj, scales)
        obj.ParameterScales = scales;
    end
    
    function fun = getCostFunction(obj)
        fun = obj.CostFunction;
    end
    
    function setCostFunction(obj, fun)
        % Set up the cost function.
        %
        % Usage
        % OPTIM.setCostFunction(FUN);
        % setCostFunction(OPTIM, FUN);
        % 
        % The input FUN can be either a function handle, or an instance of
        % CostFunction.
        
        if isa(fun, 'function_handle')
            obj.CostFunction = fun;
        elseif isa(fun, 'CostFunction')
            obj.CostFunction = @fun.evaluate;
        end
    end
    
    function fun = getOutputFunction(obj)
        fun = obj.OutputFunction;
    end
    
    function setOutputFunction(obj, fun)
        obj.OutputFunction = fun;
    end
    
    function mode = getDisplayMode(obj)
        mode = obj.DisplayMode;
    end
    
    function setDisplayMode(obj, mode)
        obj.DisplayMode = mode;
    end
    
end % general methods

%% Listeners management
methods
    function addOptimizationListener(obj, listener)
        %Adds an OptimizationListener to obj optimizer
        %
        % usage: 
        %   addOptimizationListener(OPTIM, LISTENER);
        %   OPTIM is an instance of Optimizer
        %   LISTENER is an instance of OptimizationListener
        %   The listener will listen the events of type:
        %    OptimizationStarted, 
        %    OptimizationIterated,
        %    OptimizationTerminated 
        
        % Check class of input
        if ~isa(listener, 'OptimizationListener')
            error('Input argument should be an instance of OptimizationListener');
        end
        
        % link function handles to events
        addlistener(obj, 'OptimizationStarted', @listener.optimizationStarted);
        addlistener(obj, 'OptimizationIterated', @listener.optimizationIterated);
        addlistener(obj, 'OptimizationTerminated', @listener.optimizationTerminated);
        
    end
end

end % classdef
