classdef MidRegistrationEvaluator
%MIDREGISTRATIONEVALUATOR temporary class for some test cases.
%
%   output = MidRegistrationEvaluator(input)
%
%   Example
%   MidRegistrationEvaluator
%
%   See also
%

% ------
% Author: David Legland
% e-mail: david.legland@inra.fr
% Created: 2010-11-03,    using Matlab 7.9.0.529 (R2009b)
% Copyright 2010 INRA - Cepia Software Platform.

properties
    Transfo1;
    Transfo2;
    Metric;
end

methods
    function obj = MidRegistrationEvaluator(varargin)
        
        if nargin < 3
            error('Need at least 3 arguments');
        end
        
        obj.Transfo1   = varargin{1};
        obj.Transfo2   = varargin{2};
        obj.Metric     = varargin{3};
        
    end % constructor
end

methods
    function value = evaluate(obj, params)
        
        % create parameter array for each transform
        np = length(params);
        params1 = params(1:np/2);
        params2 = params(np/2+1:end);
        
        % update each transform
        setParameters(obj.Transfo1, params1);
        setParameters(obj.Transfo2, params2);
        
        % compute resulting value
        value = computeValue(obj.Metric);
    end
end

end
