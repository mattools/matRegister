classdef WeightedSumOfFunctions < BaseFunction
%MULTIMETRICMANAGER A function that computes the sum of several functions
%
%   output = WeightedSumOfFunctions(input)
%
%   Example
%   WeightedSumOfFunctions
%
%   See also
%

% ------
% Author: David Legland
% e-mail: david.legland@inra.fr
% Created: 2010-10-27,    using Matlab 7.9.0.529 (R2009b)
% Copyright 2010 INRA - Cepia Software Platform.

%% Properties
properties
    % the set of functions
    FunctionSet;
    
    % weights associated to each function
    Weights;
    
end
 
%% Constructor
methods
    function obj = WeightedSumOfFunctions(varargin)
        
        % need 2 inputs
        if nargin<2
            error('Need 2 input arguments');
        end
        
        % check that each input is itself a ParametricFunction
        var = varargin{1};
        if ~iscell(var)
            error('Input must be a cell array of BaseFunction');
        end
        for i = 1:length(var)
            metric = var{i};
            if ~isa(metric, 'BaseFunction')
                error('WeightedSumOfFunctions:WrongClass', ...
                    'The element number %d is not a BaseFunction, but a %s', ...
                    i, class(metric));
            end
        end
        
        % stores the set of metrics
        obj.FunctionSet = var;
        
        if nargin<2
            obj.Weights = ones(1, length(obj.FunctionSet));
        else
            obj.Weights = varargin{2};
        end
        
    end % constructor
 
end % construction function
 
%% General methods
methods
 
    function value = computeValue(obj)
        % Compute the sum of the values returned by each child function
        
        % initialize
        value = 0;
        
        % compute the sum
        for i=1:length(obj.FunctionSet)
            fun = obj.FunctionSet{i};
            w = obj.Weights(i);
            
            value = value + w * computeValue(fun);
        end        
    end
    
end % general methods
 
end % classdef

