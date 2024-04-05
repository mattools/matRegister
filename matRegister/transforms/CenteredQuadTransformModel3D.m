classdef CenteredQuadTransformModel3D < CenteredTransformAbstract & ParametricTransform
% Polynomial 3D transform up to degree 2 (3*10=30 parameters).
%
%   output = CenteredQuadTransformModel3D(input)
%
%   Example
%   CenteredQuadTransformModel3D
%
%   See also
%     CenteredQuadTransformModel2D

% ------
% Author: David Legland
% e-mail: david.legland@inra.fr
% Created: 2010-11-18,    using Matlab 7.9.0.529 (R2009b)
% Copyright 2010 INRA - Cepia Software Platform.


%% Constructors
methods
    function obj = CenteredQuadTransformModel3D(varargin)
        % Create a new centered affine transform model.
        
        obj.Params = [ ...
            0 0 0 ...   % constant terms for x', y', and z'
            1 0 0 ...   % x coef 
            0 1 0 ...   % y coef 
            0 0 1 ...   % z coef 
            0 0 0 ...   % x^2 coef
            0 0 0 ...   % y^2 coef
            0 0 0 ...   % z^2 coef
            0 0 0 ...   % x*y coef
            0 0 0 ...   % x*z coef
            0 0 0 ...   % y*z coef
            ];
        
        if ~isempty(varargin)
            % extract first argument, and try to interpret
            var = varargin{1};
            if isa(var, 'CenteredQuadTransformModel3D')
                % copy constructor
                obj.Params = var.Params;
                
            elseif isnumeric(var)
                len = length(var);
                if len == 30
                    % all parameters are specified as a row vector
                    obj.Params = var;
                    
%                 elseif sum(size(var)~=[4 4])==0 || sum(size(var)~=[3 4])==0
%                     % the parameter matrix is specified
%                     var = var(1:3, :)';
%                     obj.Params = var(:)';
                    
                elseif len == 1
                    % give the possibility to call constructor with the
                    % dimension of image
                    % Initialize with identity transform
                    if var ~= 3
                        error('Defined only for 3 dimensions');
                    end
                    
                else
                    error('Please specify affine parameters');
                end
                
            else
                error('Unable to understand input arguments');
            end
        end
        
        % eventually parse additional arguments
        if nargin > 2
            if strcmp(varargin{2}, 'Center')
                % setup rotation center
                obj.Center = varargin{3};
            end
        end
        
        % setup parameter names
        obj.ParamNames = {...
            'x_w2', 'y_w2','z_w2', ...  % constant terms for x', y', and z'
            'x_wx', 'y_wx','z_wx', ...  % x coef
            'x_wy', 'y_wy','z_wy', ...  % y coef
            'x_wz', 'y_wz','z_wz', ...  % z coef
            'x_x2', 'y_x2','z_x2', ...  % x^2 coef
            'x_y2', 'y_y2','z_y2', ...  % y^2 coef
            'x_z2', 'y_z2','z_z2', ...  % z^2 coef
            'x_xy', 'y_xy','z_xy', ...  % x*y coef
            'x_xz', 'y_xz','z_xz', ...  % x*z coef
            'x_yz', 'y_yz','z_yz', ...  % x*z coef
            };
    end     
    
end


%% Methods specific to yhis class
methods    
    function initFromAffineTransform(obj, transform)
        % Initialize parameters from an affine transform (class or matrix).
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
        obj.Params = [...
            transform(1:3, 4)' ...
            transform(1:3, 1)' ...
            transform(1:3, 2)' ...
            transform(1:3, 3)' ...
            zeros(1, 18) ...
            ];
    end

end


%% Implementation of ParametricTransform methods
methods
    function jacobian = parametricJacobian(obj, x, varargin)
        % Compute jacobian matrix, i.e. derivatives for each parameter.
        
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
        x = x - obj.Center(1);
        y = y - obj.Center(2);
        z = z - obj.Center(3);
        
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
    
    function transfo = clone(obj)
        transfo = CenteredQuadTransformModel3D(obj.Params, 'Center', obj.Center);
    end
end


%% Implementation of methods inherited from Transform
methods
    function dim = getDimension(obj) %#ok<MANU>
        dim = 3;
    end

    function point2 = transformPoint(obj, point)
        % TRANSFORMPOINT Computes coordinates of transformed point.
        % PT2 = obj.transformPoint(PT);
        
        % compute centered coords.
        x = point(:, 1) - obj.Center(1);
        y = point(:, 2) - obj.Center(2);
        z = point(:, 3) - obj.Center(3);
        
        % init with translation part
        x2 = ones(size(x)) * obj.Params(1);
        y2 = ones(size(x)) * obj.Params(2);
        z2 = ones(size(x)) * obj.Params(3);
        
        % add linear contributions
        x2 = x2 + x * obj.Params(4);
        y2 = y2 + x * obj.Params(5);
        z2 = z2 + x * obj.Params(6);
        x2 = x2 + y * obj.Params(7);
        y2 = y2 + y * obj.Params(8);
        z2 = z2 + y * obj.Params(9);
        x2 = x2 + z * obj.Params(10);
        y2 = y2 + z * obj.Params(11);
        z2 = z2 + z * obj.Params(12);

        % add quadratic contributions
        x2 = x2 + x.^2 * obj.Params(13);
        y2 = y2 + x.^2 * obj.Params(14);
        z2 = z2 + x.^2 * obj.Params(15);
        x2 = x2 + y.^2 * obj.Params(16);
        y2 = y2 + y.^2 * obj.Params(17);
        z2 = z2 + y.^2 * obj.Params(18);
        x2 = x2 + z.^2 * obj.Params(19);
        y2 = y2 + z.^2 * obj.Params(20);
        z2 = z2 + z.^2 * obj.Params(21);

        % add product contributions
        x2 = x2 + x.*y * obj.Params(22);
        y2 = y2 + x.*y * obj.Params(23);
        z2 = z2 + x.*y * obj.Params(24);
        x2 = x2 + x.*z * obj.Params(25);
        y2 = y2 + x.*z * obj.Params(26);
        z2 = z2 + x.*z * obj.Params(27);
        x2 = x2 + y.*z * obj.Params(28);
        y2 = y2 + y.*z * obj.Params(29);
        z2 = z2 + y.*z * obj.Params(30);

        % recenter points
        x2 = x2 + obj.Center(1);
        y2 = y2 + obj.Center(2);
        z2 = z2 + obj.Center(3);

        % concatenate coordinates
        point2 = [x2 y2 z2];
    end
    
    function vect2 = transformVector(obj, vector, position) %#ok<INUSD,STOUT>
        % TRANSFORMVECTOR Computes coordinates of transformed vector.
        % VEC2 = obj.transformPoint(VEC, PT);
        % TODO: to be done later
        error('Not yet implemented');
    end
    
    function jacobian = jacobianMatrix(obj, point)
        % Computes jacobian matrix, i.e. derivatives wrt to each coordinate
        % jacob(i,j) = d x_i / d x_j
        
        % compute centered coords.
        x = point(:, 1) - obj.Center(1);
        y = point(:, 2) - obj.Center(2);
        z = point(:, 3) - obj.Center(3);
 
        p = obj.Params;
        dxx = p(4)  + 2*x*p(13) + y*p(22) + z*p(25);
        dyx = p(5)  + 2*x*p(14) + y*p(23) + z*p(26);
        dzx = p(6)  + 2*x*p(15) + y*p(24) + z*p(27);
        
        dxy = p(7)  + 2*y*p(16) + x*p(22) + z*p(28);
        dyy = p(8)  + 2*y*p(17) + x*p(23) + z*p(29);
        dzy = p(9)  + 2*y*p(18) + x*p(24) + z*p(30);
        
        dxz = p(10) + 2*y*p(19) + x*p(25) + y*p(28);
        dyz = p(11) + 2*y*p(20) + x*p(26) + y*p(29);
        dzz = p(12) + 2*y*p(21) + x*p(27) + y*p(30);
        
        jacobian = [dxx dxy dxz ; dyx dyy dyz ; dzx dzy dzz];
    end
    
end % Transform methods 


%% Serialization methods
methods
    function str = toStruct(obj)
        % Converts to a structure to facilitate serialization.
        str = struct('Type', 'CenteredQuadTransformModel3D', ...
            'Center', obj.Center, ...
            'Parameters', obj.Params);
        
    end
end
methods (Static)
    function transfo = fromStruct(str)
        % Creates a new instance from a structure.
        params = str.Parameters;
        transfo = CenteredQuadTransformModel3D(params, 'Center', str.Center);
    end
end


end
