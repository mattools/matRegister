% OPTIMIZERS
%
% Files
% Optimizer classes
%   Optimizer                          - General interface for single-valued function optimizers.
%   NelderMeadSimplexOptimizer         - Simplex optimizer adapted from book "Numerical Recipes 3".
%   GradientDescentOptimizer           - Gradient descent optimizer.
%   MatlabFminuncWrapper               - Encapsulation of Matlab function fminunc.
%   MatlabSimplexOptimizer             - Encapsulation of Matlab function fminsearch.
%   MultiLinearSearchOptimizer         - Optimize along successive directions.
%   BoundedMultiLinearOptimizer        - Optimize each parameter individually.
%   GaussianLinearSearchOptimizer      - Optimize each parameter individually.
%
% Functions
%   BaseFunction                       - Parent class for objects that can compute a scalar value.
%   ConstantCostFunction               - Dummy cost function used for tests.
%   CostFunction                       - Evaluate a numeric value from a parameter vector.
%   SumOfCostFunctions                 - Compute the sum of several cost functions.
%   ParametricFunction                 - Abstract classe for parametric objects that compute a value from params.
%   ParametricObject                   - Base class for all parametric objects.
%   ParametricObjectsAggregator        - Concatenates several parametric objects into a single one.
%   ParameterizedFunctionEvaluator     - Update parameters of a transform and evaluate it.
%
% Optimization monitoring
%   OptimizationListener               - Base class for listening to Optimization events.
%   OptimizationIterationLogger        - Write iteration info into a log file.
%   OptimizedValueEvolutionDisplay     - Display the evolution of optimized value.
%   ParametersEvolutionDisplay         - Display the evolution of optimization parameters.
%   ParametricFunctionEvolutionDisplay - Display the evolution of a parametric function
%
% Management of iteration step size
%   DecayFunction                      - Abstract class that models decay function for optimizers.
%   ElastixDecay                       - Decay function used in Elastix software.
%   ExponentialDecay                   - Exponential decay function.
%

% (deprecated)
%   SingleValuedCostFunction           - Evaluates a numeric value from a parameter vector
%
