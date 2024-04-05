classdef QuadPolynomialTransformModel2D < ParametricTransform
% Quadratic polynomial transform model in 2D.
%
%   Class QuadPolynomialTransformModel2D
%   Creates a new 2D quadratic transform model, with 12 parameters.
%
%   Parameters are ordered with increasing degrees in x and in y:
%   x' = P1 + P3 * x + P5 * y + P7 * x^2 +  P9 * x*y + P11 * y^2
%   y' = P2 + P4 * x + P6 * y + P8 * x^2 + P10 * x*y + P12 * y^2
%
%   Example
%   % initialize identity
%   T = QuadPolynomialTransformModel2D([0 0 1 0 0 1 zeros(1,6)]);
%   % init to a non-affine transform
%   T = QuadPolynomialTransformModel2D([50 40 1 0 0 1  [1 -2 3 -1 2 -3]*1e-3]);
%
%   See also
%

% ------
% Author: David Legland
% e-mail: david.legland@inra.fr
% Created: 2018-03-23,    using Matlab 9.3.0.713579 (R2017b)
% Copyright 2018 INRA - BIA-BIBS.


%% Properties
properties
end % end properties


%% Constructor
methods
    function obj = QuadPolynomialTransformModel2D(varargin)
    % Constructor for QuadPolynomialTransformModel2D class.
        
        % default values
        obj.Params = [ ...
            0 0 ...   % constant terms for x', y'
            1 0 ...   % x coef
            0 1 ...   % y coef
            0 0 ...   % x^2 coef
            0 0 ...   % x*y coef
            0 0 ...   % y^2 coef
            ];
        
        % initialize if necessary
        if ~isempty(varargin)
            obj.Params = varargin{1};
        end

        obj.ParamNames = {...
            'x.ww', 'y.ww', ...
            'x.xw', 'y.xw', ...
            'x.yw', 'y.yw', ...
            'x.xx', 'y.xx', ...
            'x.xy', 'y.xy', ...
            'x.yy', 'y.yy', ...
            };
    end

end % end constructors


%% Methods specific to the class
methods
    function deriv = secondDerivatives(obj, indI, indJ, point)
        % Compute the value of second derivatives for the given point(s).
        %
        % D2 = secondDerivatives(T, INDI, INDJ, POINT)
        % Return a M-by-2 array, with as many rows as the number of points.
        % First columns is the second derivative of the x-transform part,
        % and second column is the second dervative of the y-transform
        % part.
        
        % for quadratic transforms, second derivatives are constants
        if indI == 1 && indJ == 1
            deriv = obj.Params([7 8]) * 2;
        elseif indI == 2 && indJ == 2
            deriv = obj.Params([11 12]) * 2;
        elseif (indI == 1 && indJ == 2) || (indI == 2 && indJ == 1)
            deriv = obj.Params([9 10]);
        else
            error('Derivative indices must be 1 or 2');
        end
        
        % normalize size of output with number of points
        if size(point, 1) > 1
            deriv = bsxfun(@times, deriv, ones(size(point,1), 1));
        end
    end
end


%% Implementation of ParametricTransform methods
methods
    function jac = parametricJacobian(obj, x, varargin) %#ok<INUSL>
        % Compute jacobian matrix, i.e. derivatives for each parameter.
        %
        % jac = parametricJacobian(transfo, pos);
        % Returns a Jacobian matrix with 2 rows and 12 columns.
        %
        jac = [...
               1    0; ... 
               0    1; ...
               x    0; ...
               0    x; ...
               y    0; ...
               0    y; ...
            x.*x    0; ...
               0 x.*x; ...
            x.*y    0; ...
               0 x.*y; ...
            y.*y    0; ...
               0 y.*y; ...
            ]';
            
    end
    
    function transfo = clone(obj)
        transfo = QuadPolynomialTransformModel2D(obj.Params);
    end
    
end

%% Methods implementing the Transform interface
methods
    function pointT = transformPoint(obj, point)
        
        x = point(:, 1);
        y = point(:, 2);
        
        % init with translation part
        x2 = ones(size(x)) * obj.Params(1);
        y2 = ones(size(x)) * obj.Params(2);
        
        % add linear contributions
        x2 = x2 + x * obj.Params(3);
        y2 = y2 + x * obj.Params(4);
        x2 = x2 + y * obj.Params(5);
        y2 = y2 + y * obj.Params(6);
        
        % add quadratic contributions
        x2 = x2 + x.^2 * obj.Params(7);
        y2 = y2 + x.^2 * obj.Params(8);
        x2 = x2 + x.*y * obj.Params(9);
        y2 = y2 + x.*y * obj.Params(10);
        x2 = x2 + y.^2 * obj.Params(11);
        y2 = y2 + y.^2 * obj.Params(12);
        
        % concatenate coordinates
        pointT = [x2 y2];
    end
    
    function jacMat = jacobianMatrix(obj, point)
        % Computes jacobian matrix, matrix of derivatives.
        % jacob(i,j) = d x_i / d x_j
        
        % compute centered coords.
        x = point(:, 1) ;
        y = point(:, 2) ;
        
        p = obj.Params;
        dxx = p(3)  + 2*x*p(7) + y*p(9);
        dyx = p(4)  + 2*x*p(8) + y*p(10);
        
        dxy = p(5)  + 2*y*p(11) + x*p(9);
        dyy = p(6)  + 2*y*p(12) + x*p(10);
        
        jacMat = [dxx dxy ; dyx dyy];
    end
    
%     function transformVector (obj, x, varargin)
%         % Compute jacobian matrix, i.e. derivatives for each parameter
%     end
    
    function dim = getDimension (obj) %#ok<MANU>
        % Compute jacobian matrix, i.e. derivatives for each parameter.
        dim = 2;
    end
    
end % end methods


%% Serialization methods
methods
    function str = toStruct(obj)
        % Converts to a structure to facilitate serialization.
        str = struct('Type', 'QuadPolynomialTransformModel2D', ...
            'Parameters', obj.Params);
        
    end
end
methods (Static)
    function transfo = fromStruct(str)
        % Creates a new instance from a structure.
        params = str.Parameters;
        transfo = QuadPolynomialTransformModel2D(params);
    end
end

end % end classdef

