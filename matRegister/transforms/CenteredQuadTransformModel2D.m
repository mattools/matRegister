classdef CenteredQuadTransformModel2D < CenteredTransformAbstract & ParametricTransform
% Polynomial 2D transform up to degree 2 (2*6=12 parameters)
%
%   TRANSFO = CenteredQuadTransformModel2D(PARAMS)
%   Creates a new 2D centered quadratic transform model.
%
%   x' = P1 + P3 * x + P5 * y + P7 * x^2 +  P9 * x*y + P11 * y^2
%   y' = P2 + P4 * x + P6 * y + P8 * x^2 + P10 * x*y + P12 * y^2
%
%   Example
%     T = CenteredQuadTransformModel2D([0 0  1 0  0 1  .1 .1 .1 .1 .1 .1]);
%     [x y] = meshgrid(-1.5:.1:1.5, -1.5:.1:1.5);
%     pts = [x(:) y(:)];
%     figure;
%     drawPoint(pts, '.'); hold on
%     axis([-2 2 -2 2]); axis equal
%     pts2 = T.transformPoint(pts);
%     drawPoint(pts2, 'g.')
%
%   See also
%     CenteredQuadTransformModel3D
%
% ------
% Author: David Legland
% e-mail: david.legland@grignon.inra.fr
% Created: 2013-03-18,    using Matlab 7.9.0.529 (R2009b)
% Copyright 2010 INRA - Cepia Software Platform.


%% Constructors
methods
    function this = CenteredQuadTransformModel2D(varargin)
        % Create a new centered affine transform model
        
        this.params = [ ...
            0 0 ...   % constant terms for x', y', and z'
            1 0 ...   % x coef 
            0 1 ...   % y coef 
            0 0 ...   % x^2 coef
            0 0 ...   % x*y coef
            0 0 ...   % y^2 coef
            ];
        
        if ~isempty(varargin)
            % extract first argument, and try to interpret
            var = varargin{1};
            if isa(var, 'CenteredQuadTransformModel3D')
                % copy constructor
                this.params = var.params;
                
            elseif isnumeric(var)
                len = length(var);
                if len == 12
                    % all parameters are specified as a row vector
                    this.params = var;
                    
%                 elseif sum(size(var)~=[4 4])==0 || sum(size(var)~=[3 4])==0
%                     % the parameter matrix is specified
%                     var = var(1:3, :)';
%                     this.params = var(:)';
                    
                elseif len == 1
                    % give the possibility to call constructor with the
                    % dimension of image
                    % Initialize with identity transform
                    if var ~= 2
                        error('Defined only for 2 dimensions');
                    end
                    
                else
                    error('Please specify quadratic parameters');
                end
                
            else
                error('Unable to understand input arguments');
            end
        end
        
        % eventually parse additional arguments
        if nargin > 2
            if strcmp(varargin{2}, 'center')
                % setup rotation center
                this.center = varargin{3};
            end
        end
        
        % setup parameter names
        this.paramNames = {...
            'x_w2', 'y_w2', ...  % constant terms for x', y', and z'
            'x_wx', 'y_wx', ...  % x coef
            'x_wy', 'y_wy', ...  % y coef
            'x_x2', 'y_x2', ...  % x^2 coef
            'x_xy', 'y_xy', ...  % x*y coef
            'x_y2', 'y_y2', ...  % y^2 coef
            };
    end     
    
end

%% Methods specific to this class
methods    
    function initFromAffineTransform(this, transform)
        % Initialize parameters from an affine transform (class or matrix)
        %
        % Example
        % T = CenteredAffineTransformModel3D();
        % mat = createEulerAnglesRotation(.1, .2, .3);
        % T.initFromAffineTransform(mat);
        % T.getParameters()
        %
        
        % if the first argument is a Transform, extract its affine matrix
        if isa(transform, 'AffineTransform')
            transform = getAffineMatrix(transform);
        end
        
        % format matrix to have a row vector of 12 elements
        this.params = [...
            transform(1:3, 4)' ...
            transform(1:3, 1)' ...
            transform(1:3, 2)' ...
            transform(1:3, 3)' ...
            zeros(1, 18) ...
            ];
    end

end


%% Implementation of methods inherited from Transform
methods
    function dim = getDimension(this) %#ok<MANU>
        dim = 2;
    end

    function point2 = transformPoint(this, point)
        % TRANSFORMPOINT Computes coordinates of transformed point
        % PT2 = this.transformPoint(PT);
        
        % compute centered coords.
        x = point(:, 1) - this.center(1);
        y = point(:, 2) - this.center(2);
        
        % init with translation part
        x2 = ones(size(x)) * this.params(1);
        y2 = ones(size(x)) * this.params(2);
        
        % add linear contributions
        x2 = x2 + x * this.params(3);
        y2 = y2 + x * this.params(4);
        x2 = x2 + y * this.params(5);
        y2 = y2 + y * this.params(6);

        % add quadratic contributions
        x2 = x2 + x.^2 * this.params(7);
        y2 = y2 + x.^2 * this.params(8);
        x2 = x2 + x.*y * this.params(9);
        y2 = y2 + x.*y * this.params(10);
        x2 = x2 + y.^2 * this.params(11);
        y2 = y2 + y.^2 * this.params(12);

        % recenter points
        x2 = x2 + this.center(1);
        y2 = y2 + this.center(2);

        % concatenate coordinates
        point2 = [x2 y2];
    end
    
    function vect2 = transformVector(this, vector, position)
        % TRANSFORMVECTOR Computes coordinates of transformed vector
        % VEC2 = this.transformPoint(VEC, PT);
        % TODO: to be done later
        error('Not yet implemented');
    end
    
    function jacobian = getJacobian(this, point)
        % Computes jacobian matrix, i.e. derivatives wrt to each coordinate
        % jacob(i,j) = d x_i / d x_j
        
        % compute centered coords.
        x = point(:, 1) - this.center(1);
        y = point(:, 2) - this.center(2);
       
        p = this.params;
        dxx = p(3)  + 2*x*p(7) + y*p(9);
        dyx = p(4)  + 2*x*p(8) + y*p(10);
        
        dxy = p(5)  + 2*y*p(11) + x*p(9);
        dyy = p(6)  + 2*y*p(12) + x*p(10);
        
        jacobian = [dxx dxy ; dyx dyy];
    end
    
end % Transform methods 

%% Implementation of methods inherited from ParametricTransform
methods
    function jacobian = getParametricJacobian(this, x, varargin)
        % Compute jacobian matrix, i.e. derivatives for each parameter
        
        % extract coordinate of input point(s)
        if isempty(varargin)
            y = x(:,2);
            z = x(:,3);
            x = x(:,1);
        else
            y = varargin{1};
            z = varargin{2};
        end
        
        % jacobians are computed with respect to transformation center
        x = x - this.center(1);
        y = y - this.center(2);
        z = z - this.center(3);
        
        % compute the Jacobian matrix using pre-computed elements
        jacobian = [...
            1 1 1 ...
            x x x ...
            y y y ...
            z z z ...
            x.^2 x.^2 x.^2 ...
            y.^2 y.^2 y.^2 ...
            z.^2 z.^2 z.^2 ...
            x.*y x.*y x.*y ...
            x.*z x.*z x.*z ...
            y.*z y.*z y.*z ...
            ];
    end
    
end % parametric transform methods 

end
