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

