classdef NelderMeadSimplexOptimizer < Optimizer
% Simplex optimizer adapted from book "Numerical Recipes 3".
%
%   OPT = NelderMeadSimplexOptimizer()
%
%   Example
%     % Run the simplex optimizer on the Rosenbrock function
%     optimizer = NelderMeadSimplexOptimizer(@rosenbrock, [0 0], [.01 .01]);
%     [xOpt value] = optimizer.startOptimization();
%     xOpt
%       xOpt =
%           0.9993    0.9984
%     value
%       value =
%           5.5924e-006
%
%   See also
%     Optimizer, MatlabSimplexOptimizer
%

% ------
% Author: David Legland
% e-mail: david.legland@inra.fr
% Created: 2011-01-09,    using Matlab 7.9.0.529 (R2009b)
% Copyright 2010 INRA - Cepia Software Platform.

%% Properties
properties
    % maximum number of iterations
    NIters = 200;
    
    % tolerance on function value 
    FTol = 1e-5;
    
    % delta in each direction
    Deltas; 
    
    % the inner simplex, as a (ND+1)-by-ND array
    Simplex;
    
    % the vector of function evaluations
    Evals;
    
    % the sum of the vertex coordinates
    Psum;
    
    % number of function evaluation
    NumFunEvals;
end

%% Constructor
methods
    function obj = NelderMeadSimplexOptimizer(varargin)
        %Create a new Simplex Optimizer
        %
        % OPT = NelderMeadSimplexOptimizer();
        % Creates a new optimizer, whose inner fields will be initialized
        % later.
        %
        % OPT = NelderMeadSimplexOptimizer(FUN, PARAMS0);
        % OPT = NelderMeadSimplexOptimizer(FUN, PARAMS0, DELTA);
        % FUN is either a function handle, or an instance of CostFunction
        % PARAMS0 is the initial set of parameters
        % DELTA is the variation of parameters in each direction, given
        % either as a scalar, or as a row vector with the same size as the
        % parameter vector. If not specified, it is initialized as a row
        % vector the same size as PARAMS0, with values 1e-5.
        %
        
        obj = obj@Optimizer();
        
        if nargin==0
            return;
        end
        
        if nargin >= 2
            % setup cost function
            setCostFunction(obj, varargin{1});

            % setup initial optimal value
            params = varargin{2};
            setInitialParameters(obj, params);
            setParameters(obj, params);
        end

        % setup vector of variation amounts in each direction
        if nargin > 2
            del = varargin{3};
        else
            del = 1e-5;
        end
        
        if length(del) == 1
            obj.Deltas = del * ones(size(params));
        else
            obj.Deltas = del;
        end
        
    end % end of constructor
end

%% Optimizer methods

methods
    function [params, value, converged, output] = startOptimization(obj, varargin)
        %STARTOPTIMIZATION Run the optimizer, and return optimized parameters
        %
        %   PARAMS = startOptimization(OPTIM)
        %   PARAMS = OPTIM.startOptimization()
        %   Returns the optimized parameter set.
        %
        %   [PARAMS, VALUE] = startOptimization(OPTIM)
        %   [PARAMS, VALUE] = OPTIM.startOptimization()
        %   Returns the optimized parameter set and the best function evaluation.
        %
        %   [PARAMS, VALUE, CONVERGED] = startOptimization(OPTIM)
        %   [PARAMS, VALUE, CONVERGED] = OPTIM.startOptimization()
        %   Also returns a boolean indicating whether the algorithm converged or
        %   not.
        %
        %   [PARAMS, VALUE, CONVERGED, OUTPUT] = startOptimization(OPTIM)
        %   [PARAMS, VALUE, CONVERGED, OUTPUT] = OPTIM.startOptimization()
        %   Also returns a data structure containing information about the
        %   algorithm. See documentation of 'fminsearch' for details.
        %
        %
        %   Example
        %     startOptimization
        %
        %   See also
        %
        
        % ------
        % Author: David Legland
        % e-mail: david.legland@inra.fr
        % Created: 2011-01-09,    using Matlab 7.9.0.529 (R2009b)
        % Copyright 2010 INRA - Cepia Software Platform.
        
        %TODO: provide psb to start optimization with a specified simplex
        
        TINY = 1e-10;
        
        
        %% Initialization
        
        % setup initial parameters
        params = obj.Params;
        if ~isempty(obj.InitialParameters)
            params = obj.InitialParameters;
        end
        if ~isempty(varargin)
            params = varargin{1};
        end
        obj.Params = params;
        
        % Notify beginning of optimization
        notify(obj, 'OptimizationStarted');
        
        % initialize the simplex.
        initializeSimplex(obj);
        
        % state of the algorithm
        exitMessage = 'Algorithm started';
        converged = false;
        
        
        %% Main loop
        
        % infinite loop
        iter = 1;
        while true
            % first, determines the indices of points with the highest (i.e.
            % worst), next highest, and lowest (i.e. best) values.
            [dummy, indices] = sort(obj.Evals); %#ok<ASGLU>
            indLow  = indices(1);
            indHigh = indices(end);
            indNext = indices(end-1);
            
            % update optimized value and position
            obj.Params = obj.Simplex(indLow, :);
            obj.Value  = obj.Evals(indLow);
            
            % compute relative difference between highest and lowest
            fLow    = obj.Evals(indLow);
            fHigh   = obj.Evals(indHigh);
            rtol = 2 * abs(fHigh - fLow) / (abs(fHigh) + abs(fLow) + TINY);
            
            % termination with function evaluation
            if rtol < obj.FTol
                exitMessage = sprintf('Function converged with relative tolerance %g', obj.FTol);
                converged = true;
                break;
            end
            
            % begin a new iteration
            
            % first extrapolate by a factor -1 through the face of the simplex
            % opposite to the highest point.
            [xTry, fTry] = evaluateReflection(obj, indHigh, -1);
            
            % if the value at the evaluated position is better than current
            % highest value, then replace the highest value
            if fTry < obj.Evals(indHigh)
                updateSimplex(obj, indHigh, xTry, fTry);
                if strcmp(obj.DisplayMode, 'iter')
                    disp('reflection');
                end
                notify(obj, 'OptimizationIterated');
            end
            
            if fTry <= obj.Evals(indLow)
                % if new evaluation is better than current minimum, try to expand
                [xTry, fTry] = evaluateReflection(obj, indHigh, 2);
                
                if fTry < obj.Evals(indHigh)
                    % expansion was successful
                    updateSimplex(obj, indHigh, xTry, fTry);
                    if strcmp(obj.DisplayMode, 'iter')
                        disp('expansion');
                    end
                    notify(obj, 'OptimizationIterated');
                end
                
            elseif fTry >= obj.Evals(indNext)
                % if new evaluation is worse than the second-highest point, look
                % for an intermediate point (i.e. do a one-dimensional contraction)
                [xTry, fTry] = evaluateReflection(obj, indHigh, .5);
                
                if fTry < obj.Evals(indHigh)
                    % contraction was successful
                    updateSimplex(obj, indHigh, xTry, fTry);
                    if strcmp(obj.DisplayMode, 'iter')
                        disp('contraction');
                    end
                    notify(obj, 'OptimizationIterated');
                    
                else
                    % 1D contraction was not successful, so perform a shrink
                    % (multidimensional contraction) around lowest point
                    contractSimplex(obj, indLow);
                    if strcmp(obj.DisplayMode, 'iter')
                        disp('shrink');
                    end
                    notify(obj, 'OptimizationIterated');
                end
                
            end
            
            % termination with number of iterations
            if iter > obj.NIters
                exitMessage = sprintf('Iteration number reached maximum allowed value: %d', obj.nIter);
                break;
            end
            
            iter = iter + 1;
        end % main iteration loop
        
        
        %% Terminates
        
        if converged
            if strcmp(obj.DisplayMode, {'iter', 'final'})
                disp(exitMessage);
            end
        else
            if strcmp(obj.DisplayMode, {'iter', 'final', 'notify'})
                disp(exitMessage);
            end
        end
        
        % send termination event
        notify(obj, 'OptimizationTerminated');
        
        % return the current state of the optimizer
        params  = obj.Params;
        value   = obj.Value;
        
        % create output structure
        output.algorithm    = '';
        output.funcCount    = obj.NumFunEvals;
        output.iterations   = iter;
        output.message      = exitMessage;
        
    end
end

%% Private functions
methods (Access = private)
    function initializeSimplex(obj)
        % Initialize the simplex. 
        % This is a (ND+1)-by-ND array containing the coordinates of a
        % vertex on each row.
        
        if isempty(obj.Params)
            if isempty(obj.InitialParameters)
                error('MatRegister:NelderMeadSimplexOptimizer:initialization', ...
                    'Need to specify initial parameters');
            end
            obj.Params = obj.InitialParameters;
        end
        
        nd = length(obj.Params);
        
        % ensure delta has valid value
        if isempty(obj.Deltas)
            obj.Deltas = 1e-5 * ones(1, nd);
        end
        
        % initialize vertex coordinates
        obj.Simplex = repmat(obj.Params, nd+1, 1);
        for i = 1:nd
            obj.Simplex(i+1, i) = obj.Params(i) + obj.Deltas(i);
        end
        
        % compute sum of vertex coordinates
        obj.Psum = sum(obj.Simplex, 1);

        % evaluate function for each vertex of the simplex
        obj.Evals = zeros(nd+1, 1);
        for i=1:nd+1
            obj.Evals(i) = obj.CostFunction(obj.Simplex(i, :));
        end
        
        obj.NumFunEvals = 0;

   end
    
    function [ptry, ytry] = evaluateReflection(obj, ihi, fac)
        % Helper function that evaluates the value of the function at the
        % reflection of point with index ihi
        %
        % [PTRY YTRY] = evaluateReflection(PT_INDEX, FACTOR)
        % PT_INDEX index of the vertex that is updated
        % FACTOR expansion (>1), reflection (<0) or contraction (0<F<1)
        %   factor 
        % PTRY is the new computed point
        % YTRY is the function evaluation at the newly evaluated point
        
        % compute weighting factors
        nd = length(obj.Params);
        fac1 = (1 - fac) / nd;
        fac2 = fac1 - fac;
         
        % position of the new candidate point
        ptry = obj.Psum * fac1 - obj.Simplex(ihi, :) * fac2;
        
        % evaluate function value
        ytry = obj.CostFunction(ptry);
        obj.NumFunEvals = obj.NumFunEvals + 1;
    end

    function updateSimplex(obj, ihi, pTry, yTry)
        obj.Evals(ihi) = yTry;
        obj.Psum = obj.Psum - obj.Simplex(ihi, :) + pTry ;
        obj.Simplex(ihi, :) = pTry;
    end
    
    function contractSimplex(obj, indLow)
        nd = length(obj.Params);
        pLow = obj.Simplex(indLow,:);
        for i = [1:indLow-1 indLow+1:nd]
            obj.Simplex(i, :) = (obj.Simplex(i,:) + pLow) * .5;
            obj.Evals(i) = obj.CostFunction(obj.Simplex(i,:));
        end
        
        obj.NumFunEvals = obj.NumFunEvals + nd;
    end
    
end % private methods

end % classdef