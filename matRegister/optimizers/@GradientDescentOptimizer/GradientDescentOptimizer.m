classdef GradientDescentOptimizer < Optimizer
%GRADIENTDESCENTOPTIMIZER Gradient descent optimizer
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

end
