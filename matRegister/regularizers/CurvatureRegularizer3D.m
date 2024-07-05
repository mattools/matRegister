classdef CurvatureRegularizer3D < TransformRegularizer
%CURVATUREREGULARIZER2D Curvature regularizer for 3D transforms.
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
% e-mail: david.legland@inrae.fr
% Created: 2018-08-08,    using Matlab 9.4.0.813654 (R2018a)
% Copyright 2018 INRA - BIA-BIBS.


%% Properties
properties
    Transform;
    
end % end properties


%% Constructor
methods
    function obj = CurvatureRegularizer3D(varargin)
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
        
        d11 = secondDerivatives(obj.Transform, points, 1, 1);
        d22 = secondDerivatives(obj.Transform, points, 2, 2);
        d33 = secondDerivatives(obj.Transform, points, 3, 3);
        val = sum(d11.^2 + d22.^2 + d33.^2, 2);
    end
end % end methods

end % end classdef

