classdef AffineTransform < Transform
%AFFINETRANSFORM  Abstract class for AffineTransform.
%   
%   AffineTransform are designed to transform points into other points
%   using linear relations.
%
%   T : R^d -> R^d
%   
%   few methods are already implemented for convenience:
%   - transformPoint
%   - transformVector
%   - jacobianMatrix
%   They all use the abstract getAffineMatrix method.
%   
%
%   See also
%   Transform, MatrixAffineTransform
%

% ------
% Author: David Legland
% e-mail: david.legland@inra.fr
% Created: 2010-04-09,    using Matlab 7.9.0.529 (R2009b)
% Copyright 2010 INRA - Cepia Software Platform.

properties
    % the set of inner parameters of the AffineTransform
end

%% Static methods
methods(Static)
    function transfo = createTranslation(shift)
        % Creates a translation from a shift vector
        %
        % Example:
        % transfo = AffineTransform.createTranslation([5 3]);
        
        nd = length(shift);
        mat = [eye(nd) shift(:) ; zeros(1, nd) 1];
        transfo = MatrixAffineTransform(mat);
    end
    
    function transfo = createScaling(factors)
        % Creates a scaling transform from a set of scaling factors
        %
        % Example:
        % transfo = AffineTransform.createScaling([3 2]);
        
        nd = length(factors);
        mat = [diag(factors) zeros(nd, 1) ; zeros(1, nd) 1];
        transfo = MatrixAffineTransform(mat);
    end
    
    function transfo = createRotation(varargin)
        % Creates a 2D rotation around origin, from an angle in radians
        
        if nargin == 1
            cx = 0; cy = 0;
            theta = varargin{1};
        elseif nargin == 2
            var1 = varargin{1};
            cx = var1(1);
            cy = var1(2);
        end
        cot = cos(theta);
        sit = sin(theta);
        tx =  cy*sit - cx*cot + cx;
        ty = -cy*cot - cx*sit + cy;
        
        % create transformation matrix
        mat = [cot -sit tx; sit cot ty; 0 0 1];
        transfo = MatrixAffineTransform(mat);
    end
    
    function transfo = createRotationOx(varargin)
        % Creates a 3D rotation around X axis, from an angle in radians
        %
        % Usage:
        % transfo = AffineTransform.createRotationOx(angle);
        % transfo = AffineTransform.createRotationOx(center, angle);
        
        if nargin == 1
            cx = 0; cy = 0; cz = 0;
            theta = varargin{1};
        elseif nargin == 2
            var1 = varargin{1};
            cx = var1(1);
            cy = var1(2);
            cz = var1(3);
        end
        
        % compute rotation part
        cot = cos(theta);
        sit = sin(theta);
        mat = [...
            1 0 0 0;...
            0 cot -sit 0;...
            0 sit cot 0;...
            0 0 0 1];

        % recenter transform
        t = [1 0 0 cx;0 1 0 cy;0 0 1 cz;0 0 0 1];
        mat = t * mat / t; 
        
        % create transform instance
        transfo = MatrixAffineTransform(mat);
    end
    
    function transfo = createRotationOy(varargin)
        % Creates a 3D rotation around Y axis, from an angle in radians
        %
        % Usage:
        % transfo = AffineTransform.createRotationOy(angle);
        % transfo = AffineTransform.createRotationOy(center, angle);
        
        if nargin == 1
            cx = 0; cy = 0; cz = 0;
            theta = varargin{1};
        elseif nargin == 2
            var1 = varargin{1};
            cx = var1(1);
            cy = var1(2);
            cz = var1(3);
        end
        
        % compute rotation part
        cot = cos(theta);
        sit = sin(theta);
        mat = [...
            cot  0  sit  0;...
            0    1    0  0;...
            -sit 0  cot  0;...
            0    0    0  1];

        % recenter transform
        t = [1 0 0 cx;0 1 0 cy;0 0 1 cz;0 0 0 1];
        mat = t * mat / t; 
        
        % create transform instance
        transfo = MatrixAffineTransform(mat);
    end
    
    function transfo = createRotationOz(varargin)
        % Creates a 3D rotation around Z axis, from an angle in radians
        %
        % Usage:
        % transfo = AffineTransform.createRotationOz(angle);
        % transfo = AffineTransform.createRotationOz(center, angle);
        
        if nargin == 1
            cx = 0; cy = 0; cz = 0;
            theta = varargin{1};
        elseif nargin == 2
            var1 = varargin{1};
            cx = var1(1);
            cy = var1(2);
            cz = var1(3);
        end
        
        % compute rotation part
        cot = cos(theta);
        sit = sin(theta);
        mat = [...
            cot -sit 0 0;...
            sit  cot 0 0;...
            0 0 1 0;...
            0 0 0 1];

        % recenter transform
        t = [1 0 0 cx;0 1 0 cy;0 0 1 cz;0 0 0 1];
        mat = t * mat / t; 
        
        % create transform instance
        transfo = MatrixAffineTransform(mat);
    end
    
end

%% Abstract methods
methods(Abstract)
    affineMatrix(obj)
    % Returns the (d+1)*(d+1) affine matrix that represents obj transform
    %
    % usage:
    % trans = Rotation2D([3 4], pi/3); % define an affine transform
    % mat = getAffineMatrix(trans);   % extract matrix
end

%% General methods specific to Affine transforms
methods
    function mat = getAffineMatrix(obj)
        warning('deprecated: use method affineMatrix instead');
        mat = affineMatrix(obj);
    end
    
    function res = mtimes(obj, that)
        % Multiplies two affine transforms
        % RES = THIS * THAT
        % RES = mtimes(THIS*THAT)
        % both THIS and THAT must be affine transforms.
        %
        % the result is an instance of MatrixAffineTransform.
        
        % extract matrices
        mat1 = affineMatrix(obj);
        mat2 = affineMatrix(that);
        
        % check sizes are equal
        if sum(size(mat1) ~= size(mat2)) > 0
            error('The two transforms must have matrices the same size');
        end
        
        % returns the multipled transform
        res = MatrixAffineTransform(mat1 * mat2);
    end
    
    function res = getInverse(obj)
        % Computes the inverse transform of obj affine transform
        % 
        % TINV = T.getInverse();
        % or 
        % TINV = getInverse(T);
        %
        mat = affineMatrix(obj);
        res = MatrixAffineTransform(inv(mat));
    end
end

%% Methods implementing the transform interface
methods
    function dim = getDimension(obj)

        % extract affine coefficients
        mat = affineMatrix(obj);
        dim = size(mat, 1) - 1;
    end

    function varargout = transformPoint(obj, point, varargin)
        % Transform the point using affine coefficient
        % 
        % transfo.transformPoint(POINT);
        % transfo.transformPoint(X, Y);
        % transfo.transformPoint(X, Y, Z);
        
        % extract affine coefficients
        mat = affineMatrix(obj);
        
        % format to process a single array
        baseSize = size(point);
        if ~isempty(varargin)
            point = point(:);
            point(1, nargin-1) = 0;
            for i = 3:nargin
                var = varargin{i-2};
                point(:,i) = var(:);
            end
        end

        if size(mat, 2)-1 ~= size(point, 2)
            error('Point and transform must have same dimension');
        end
        
        % compute coordinate of result point
        res = zeros(size(point));
        for i = 1:size(point, 2)
            res(:,i) = point * mat(i, 1:end-1)' + mat(i, end);
        end
        
        % format output arguments
        if nargout <= 1
            varargout = {res};
        else
            varargout = cell(1, nargout);
            for i = 1:nargout
                varargout{i} = reshape(res(:,i), baseSize);
            end
        end
    end
    
    function vector = transformVector(obj, vector, position) %#ok<INUSD>
        % Transform the vector using affine coefficients
        % 
        % VEC2 = transfo.transformVector(VEC, POINT);

        mat = affineMatrix(obj);
        res = zeros(size(point));
        for i = 1:size(point, 2)
            res(:,i) = point * mat(i, 1:end-1)';
        end
    end

    function jacMat = jacobianMatrix(obj)
        % Compute jacobian matrix, i.e. derivatives for coordinate
        % jacob(i,j) = d x_i / d x_j
        mat = affineMatrix(obj);
        jacMat = mat(1:end-1, 1:end-1);
    end

end % methods

end % classdef
