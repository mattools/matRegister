classdef GradientDescentOptimizer < Optimizer
%@GRADIENTDESCENTOPTIMIZER Gradient descent optimizer
%
%   output = @GradientDescentOptimizer(input)
%
%   Example
%   @GradientDescentOptimizer
%
%   See also
%
%
% ------
% Author: David Legland
% e-mail: david.legland@grignon.inra.fr
% Created: 2010-10-07,    using Matlab 7.9.0.529 (R2009b)
% Copyright 2010 INRA - Cepia Software Platform.

properties
    %TODO: add gradient descent control
    nIter   = 200;
    
    % the decay function that controls the length of displacement
    decayFunction = ExponentialDecay(200);
    
    % current iteration value
    iter;
    
    % current value of function gradient
    gradient;
end

%% Constructor
methods
    function this = GradientDescentOptimizer(varargin)
        % Create a new Gradient descent optimizer
        %
        % Default constructor is empty constructor. It is also possible to
        % initialize with cost function and initial set of parameters.
        % 
        %   Example
        %   OPTIM = GradientDescentOptimizer();
        %   setDecayFunction(OPTIM, ExponentialDecay(50));
        %
        
        % call the parent constructor
        this = this@Optimizer(varargin{:});
    end
    
end


%% General methods
methods
    function setDecayFunction(this, decayFcn)
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
        
        this.decayFunction = decayFcn;
    end
end

end
