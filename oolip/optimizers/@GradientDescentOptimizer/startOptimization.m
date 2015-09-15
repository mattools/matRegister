function [params, value] = startOptimization(this, varargin)
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
% e-mail: david.legland@grignon.inra.fr
% Created: 2010-10-07,    using Matlab 7.9.0.529 (R2009b)
% Copyright 2010 INRA - Cepia Software Platform.


%% Initialisation

% optimisation parameters
nIter   = this.nIter;

% allocate memory
step = zeros(nIter, 1);

% setup parameters to initial value
if ~isempty(this.initialParameters)
    this.params = this.initialParameters;
end
if ~isempty(varargin)
    this.params = varargin{1};
end

% initialize optimization result
this.bestValue = inf;
this.bestParams = this.params;

% reset direction vector to null vector
% -> first iteration will be equal to the initial parameter vector
direction = zeros(size(this.params));

% Notify beginning of optimization
this.notify('OptimizationStarted');


for i = 1:nIter
    %% Update parameters and dependent parameters

    this.iter = i;
    
    % compute step depending on current iteration
    step(i) = evaluate(this.decayFunction, i);
    
    % compute new set of parameters
    this.params = this.params + direction * step(i);

    % update metric
    [this.value, this.gradient] = this.costFunction(this.params);
    
    % if value is better than before, update best value
    if this.value < this.bestValue
        this.bestValue = this.value;
        this.bestParams = this.params;
    end
    
    % if scales are initialized, scales the derivative
    if ~isempty(this.parameterScales)
        % dimension check
        if length(this.parameterScales) ~= length(this.params)
            error('Scaling parameters should have same size as parameters');
        end
        
        this.gradient = this.gradient ./ this.parameterScales;
    end
    
    % search direction (with a minus sign because we are looking for the
    % minimum)
    direction = -this.gradient / norm(this.gradient);
    
    
    %% Notifications   
    
    % Notify the end of iteration to OptimizationListeners
    this.notify('OptimizationIterated');

    % Call an output function for processing about current point
    % (for compatibility with Matlab syntax)
    if ~isempty(this.outputFunction)
        % setup optim values
        optimValues.fval = this.value;
        optimValues.iteration = i;
        optimValues.procedure = 'Gradient descent';
        
        % call output function with appropriate parameters
        stop = this.outputFunction(this.params, optimValues, 'iter');
        if stop
            break;
        end
    end
end


%% Finalisation

% Notify termination event
this.notify('OptimizationTerminated');

% returns the current set of parameters
value = this.bestValue;
params = this.bestParams;

