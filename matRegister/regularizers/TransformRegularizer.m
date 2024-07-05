classdef TransformRegularizer < handle
%TRANSFORMREGULARIZER Evaluate the regularisation function of a transform
%
%   Class TransformRegularizer
%
%   Example
%   TransformRegularizer
%
%   See also
%

% ------
% Author: David Legland
% e-mail: david.legland@inra.fr
% Created: 2018-08-08,    using Matlab 9.4.0.813654 (R2018a)
% Copyright 2018 INRA - BIA-BIBS.


%% Properties
properties
end % end properties


%% Constructor
methods
    function obj = TransformRegularizer(varargin)
    % Constructor for TransformRegularizer class

    end

end % end constructors


%% Methods
methods (Abstract)
    % Evaluate the regularization for a given point.
    evaluate(obj, point);
end % end methods

end % end classdef

