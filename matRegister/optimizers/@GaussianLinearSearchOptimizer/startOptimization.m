function [params, value] = startOptimization(obj, varargin)
%STARTOPTIMIZATION  Run the optimizer, and return optimized parameters
%
%   PARAMS = startOptimization(OPTIM)
%   PARAMS = OPTIM.startOptimization()
%
%   [PARAMS VALUE] = startOptimization(OPTIM)
%   [PARAMS VALUE] = OPTIM.startOptimization()
%
%   Example
%   StartOptimization
%
%   See also
%

% ------
% Author: David Legland
% e-mail: david.legland@inra.fr
% Created: 2010-11-24,    using Matlab 7.9.0.529 (R2009b)
% Copyright 2010 INRA - Cepia Software Platform.

% Notify beginning of optimization
notify(obj, 'OptimizationStarted');

% get initial parameters
params = obj.Params;
if ~isempty(obj.InitialParameters)
    params = obj.InitialParameters;
end
if ~isempty(varargin)
    params = varargin{1};
end

% extract variability
nParams = length(params);
variab  = obj.ParameterVariability;

if isscalar(variab)
    variab = repmat(variab, 1, nParams);
end
if length(variab) ~= nParams
    warning('oolip:NonInitializedParameter',...
        'Variability not specified, use default variability');
    variab = ones(size(params));
    obj.ParameterVariability = variab;
end

% generate nValues values between 0 and 1, avoiding bounds
pv = linspace(0, 1, obj.NValues+2);
pv = pv(2:end-1);

for i = 1:obj.NIters
    if strcmp(obj.DisplayMode, 'iter')
        disp(sprintf('iter %d/%d', i, obj.NIters)); %#ok<DSPS>
    end

    for p = 1:nParams
        % choose a set of values distributed around the current value
        % with a gaussian distribution
        paramValues = norminv(pv, params(p), variab(p));
        
        % compute the metric for the given param
        res = zeros(obj.NValues, 1);
        for k = 1:obj.NValues
            params(p) = paramValues(k);
            res(k) = obj.CostFunction(params);
        end
        
        [value, bestK] = min(res);
        params(p) = paramValues(bestK);
        
        % update optimizer internal state
        obj.Params = params;
        obj.Value  = value;
        
        % notify
        notify(obj, 'OptimizationIterated');
    end 
end

% update inner data
obj.Params = params;
obj.Value = value;

% Notify the end of optimization
notify(obj, 'OptimizationTerminated');

