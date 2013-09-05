% OPTIMIZERS
%
% Files
% Optimizers
%   Optimizer                          - General interface for single-valued function optimizers
%   BoundedMultiLinearOptimizer        - Optimize each parameter individually
%   GaussianLinearSearchOptimizer      - Optimize each parameter individually
%   GradientDescentOptimizer           - Gradient descent optimizer
%   MatlabFminuncWrapper               - One-line description here, please.
%   MatlabSimplexOptimizer             - Encapsulation of Matlab function fminsearch
%   MultiLinearSearchOptimizer         - Optimize along successive directions
%   NelderMeadSimplexOptimizer         - Simplex optimizer adapted from Numerical Recipes 3
%
% Functions
%   BaseFunction                       - Parent class for objects that can compute a scalar value
%   ConstantCostFunction               - Dummy cost function used for tests
%   CostFunction                       - Evaluate a numeric value from a parameter vector
%   SingleValuedCostFunction           - Evaluates a numeric value from a parameter vector
%   SumOfCostFunctions                 - Compute the sum of several cost functions
%   ParametricFunction                 - One-line description here, please.
%   ParametricObject                   - Base class for all parametric objects
%   ParametricObjectsAggregator        - Concatenates several parametric objects into one
%   ParameterizedFunctionEvaluator     - One-line description here, please.
%
% Optimizaion utilities
%   DecayFunction                      - Abstract class that models decay function for optimizers
%   ElastixDecay                       - Decay function used in Elastix software
%   ExponentialDecay                   - Exponential decay function.
%   OptimizationIterationLogger        - One-line description here, please.
%   OptimizationListener               - Base class for listening to Optimization events
%   OptimizedValueEvolutionDisplay     - Displays evolution of optimized value
%   ParametersEvolutionDisplay         - Displays evolution of optmization parameters
%   ParametricFunctionEvolutionDisplay - Displays evolution of a parametric function
