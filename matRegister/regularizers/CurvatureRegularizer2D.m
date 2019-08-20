classdef CurvatureRegularizer2D < TransformRegularizer
%CURVATUREREGULARIZER2D Curvature regularizer for 2D transforms
%
%   Class CurvatureRegularizer2D
%
%   Example
%     T = QuadPolynomialTransformModel2D([50 40 1 0 0 1  [1 -2 3 -1 2 -3]*1e-3]);
%     lx = -50:5:50;
%     [x, y] = meshgrid(lx, lx);
%     pts = [x(:) y(:)];
%     regul = CurvatureRegularizer2D(T);
%     values = evaluate(regul, pts);
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
    Transform;
    
end % end properties


%% Constructor
methods
    function obj = CurvatureRegularizer2D(varargin)
    % Constructor for CurvatureRegularizer2D class
    
        if ~isa(varargin{1}, 'Transform')
            error('Requires Transform class as first argument');
        end
        obj.Transform = varargin{1};

    end

end % end constructors


%% Methods implementing the TransformRegularizer interface
methods
    function val = evaluate(obj, points)
        
        d11 = secondDerivatives(obj.Transform, 1, 1, points);
        d22 = secondDerivatives(obj.Transform, 2, 2, points);
        val = sum(d11.^2 + d22.^2, 2);
    end
end % end methods

end % end classdef

