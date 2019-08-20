function [params, value] = startOptimization(obj)
%STARTOPTIMIZATION Run the optimization algorithm


% ensure there is a valid direction set
if isempty(obj.DirectionSet)
    initDirectionSet(obj);
end


% default tolerance
tol = 1e-5;

nDirs = size(obj.DirectionSet, 1);

% main loop
% Use linear search in a set of directions that span the parameters of the
% first transform.
% The same process could be obteained by calling 'directionSetMinimizer'.
dirIndex = 1;
for i = 1:obj.NIters
    if obj.Verbose
        disp(sprintf('iteration %d / %d', i, obj.NIters)); %#ok<DSPS>
    end
    
    % current direction
    dir = obj.DirectionSet(dirIndex, :);
    
    % update direction index
    dirIndex = mod(dirIndex, nDirs) + 1;
    
    % use a function handle of 1 variable
    fun1 = @(t) obj.CostFunction(obj.Params + t*dir);
    
    % guess initial bounds of the function
    ax = 0;
    bx = 1;
    [ax, bx, cx] = fMinBracket(fun1, ax, bx);
    
    % search minimum along dimension DIR
    [tmin, value] = brentLineSearch(fun1, ax, bx, cx, tol);
    
    % construct new optimal point
    obj.Params = obj.Params + tmin*dir;
    
    % Call an output function for processing current point
    if ~isempty(obj.OutputFunction)
        % setup optim values
        optimValues.fval = value;
        optimValues.iteration = i;
        optimValues.procedure = 'Linear search';
        
        % call output function with appropriate parameters
        stop = obj.OutputFunction(obj.Params, optimValues, 'iter');
        if stop
            break;
        end
    end
end

% format output arguments
params = obj.Params;

% we need a final end because nested functions are used
end
