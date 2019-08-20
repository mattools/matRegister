classdef GradientDescentOptimizer < Optimizer
% Gradient descent optimizer.
%
%   OPTIM = GradientDescentOptimizer();
%   Initializes an empty optimizer.
%
%   OPTIM = GradientDescentOptimizer(FUN, X0);
%   Initializes gradient desczent optimizer with a cost function FUN (given
%   as a function handle) and a starting point X0 (given as a 1-by-P row
%   vector).
%
%   Example
%     x0 = [0 0];
%     optim = GradientDescentOptimizer(@rosenbrock, x0);
%     setDecayFunction(optim, ExponentialDecay(200));
%     optim.nIter = 1000;
%     xHat = optim.startOptimization([0 0])
%     xHat =
%         [1.0000 0.9999]
%
%
%   See Also
%     Optimizer, rosenbrock
%

% ------
% Author: David Legland
% e-mail: david.legland@inra.fr
% Created: 2010-10-07,    using Matlab 7.9.0.529 (R2009b)
% Copyright 2010 INRA - Cepia Software Platform.

properties
    %TODO: add gradient descent control
    NIters   = 200;
    
    % the decay function that controls the length of displacement
    DecayFunction = ExponentialDecay(200);
    
    % current iteration value
    Iter;
    
    % current value of function gradient
    Gradient;
    
    BestValue;
    BestParams;
end

%% Constructor
methods
    function obj = GradientDescentOptimizer(varargin)
        % Create a new Gradient descent optimizer
        %
        % Default constructor is empty constructor. It is also possible to
        % initialize with cost function and initial set of parameters.
        % 
        %   Example
        %   OPTIM = GradientDescentOptimizer();
        %   setDecayFunction(OPTIM, ExponentialDecay(50));
        %
        %   fun = @(x) (x(1)-4)^2 + (x(2)-3)^2;
        %   x0 = [1 1];
        %   optim = GradientDescentOptimizer(fun, x0);
        %   xHat = optim.startOptimization();
        %
        
        % call the parent constructor
        obj = obj@Optimizer(varargin{:});
    end
    
end


%% General methods
methods
    function setDecayFunction(obj, decayFcn)
        % Changes the decay function
        %
        %   setDecayFunction(OPTIM, DECAY)
        %   Replaces the current decay function of the optimizer by DECAY.
        %   DECAY must be an instance of DecayFunction.
        %
        %   Example
        %     setDecayFunction(OPTIM, ExponentialDecay(50));
        %
        
        if ~isa(decayFcn, 'DecayFunction')
            error('decay function should be a subclass of "DecayFunction"');
        end
        
        obj.DecayFunction = decayFcn;
    end
end


%% Optimization Methods

methods
    function [params, value] = startOptimization(obj, varargin)
        %STARTOPTIMIZATION Start the gradient descent optimization algorithm
        %
        %   xHat = startOptimization(OPTIM)
        %   xHat = startOptimization(OPTIM, X0)
        %
        %   Example
        %   startOptimization
        %
        %   See also
        %
        
        % ------
        % Author: David Legland
        % e-mail: david.legland@inra.fr
        % Created: 2010-10-07,    using Matlab 7.9.0.529 (R2009b)
        % Copyright 2010 INRA - Cepia Software Platform.
        
        
        %% Initialisation
        
        % optimisation parameters
        nIter = obj.NIters;
        
        % allocate memory
        step = zeros(nIter, 1);
        
        % setup parameters to initial value
        if ~isempty(obj.InitialParameters)
            obj.Params = obj.InitialParameters;
        end
        if ~isempty(varargin)
            obj.Params = varargin{1};
        end
        
        % initialize optimization result
        obj.BestValue = inf;
        obj.BestParams = obj.Params;
        
        % reset direction vector to null vector
        % -> first iteration will be equal to the initial parameter vector
        direction = zeros(size(obj.Params));
        
        % Notify beginning of optimization
        notify(obj, 'OptimizationStarted');
        
        
        for i = 1:nIter
            %% Update parameters and dependent parameters
            
            obj.Iter = i;
            
            % compute step depending on current iteration
            step(i) = evaluate(obj.DecayFunction, i);
            
            % compute new set of parameters
            obj.Params = obj.Params + direction * step(i);
            
            % update metric
            [obj.Value, obj.Gradient] = obj.CostFunction(obj.Params);
            
            % if value is better than before, update best value
            if obj.Value < obj.BestValue
                obj.BestValue = obj.Value;
                obj.BestParams = obj.Params;
            end
            
            % if scales are initialized, scales the derivative
            if ~isempty(obj.ParameterScales)
                % dimension check
                if length(obj.ParameterScales) ~= length(obj.Params)
                    error('Scaling parameters should have same size as parameters');
                end
                
                obj.Gradient = obj.Gradient ./ obj.ParameterScales;
            end
            
            % search direction (with a minus sign because we are looking for the
            % minimum)
            direction = -obj.Gradient / norm(obj.Gradient);
            
            
            %% Notifications
            
            % Notify the end of iteration to OptimizationListeners
            notify(obj, 'OptimizationIterated');
            
            % Call an output function for processing about current point
            % (for compatibility with Matlab syntax)
            if ~isempty(obj.OutputFunction)
                % setup optim values
                optimValues.fval = obj.Value;
                optimValues.Iteration = i;
                optimValues.procedure = 'Gradient descent';
                
                % call output function with appropriate parameters
                stop = obj.OutputFunction(obj.Params, optimValues, 'iter');
                if stop
                    break;
                end
            end
        end
        
        
        %% Finalisation
        
        % Notify termination event
        notify(obj, 'OptimizationTerminated');
        
        % returns the current set of parameters
        value = obj.BestValue;
        params = obj.BestParams;
        
    end
end

end
