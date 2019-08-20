classdef MutualInformationMetric < ImageToImageMetric
% Compute minus the mutual information of two images.
%
%   output = MutualInformationMetric(input)
%
%   Example
%   METRIC = MutualInformationMetric(IMG1, IMG2, POINTS);
%   METRIC.computeValue()
%
%   See also
%

% ------
% Author: David Legland
% e-mail: david.legland@inra.fr
% Created: 2010-08-12,    using Matlab 7.9.0.529 (R2009b)
% Copyright 2010 INRA - Cepia Software Platform.

%% Constructor
methods
    function obj = MutualInformationMetric(varargin)
        % calls the parent constructor
        obj = obj@ImageToImageMetric(varargin{:});
        
    end % constructor
    
end % methods

end % classdef
