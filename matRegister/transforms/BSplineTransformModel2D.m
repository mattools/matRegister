classdef BSplineTransformModel2D < ParametricTransform
%BSPLINETRANSFORMMODEL2D Cubic B-Spline Transform model in 2D.
%
%   Class BSplineTransformModel2D
%
%   Grid is composed of M-by-N vertices, with M number of rows and N number
%   of columns. Iteration along x direction first.
%   Parameters correspond to shift vector associated to each vertex:
%   [vx11 vy11 vx12 vy12 ... vxIJ vyIJ ... vxMN vyMN]
%
%   Example
%   BSplineTransformModel2D
%
%   See also
%     BSplineTransformModel3D

% ------
% Author: David Legland
% e-mail: david.legland@inra.fr
% Created: 2018-08-09,    using Matlab 9.4.0.813654 (R2018a)
% Copyright 2018 INRA - BIA-BIBS.


%% Properties
properties
    % Number of vertices of the grid in each direction
    % (as a 1-by-2 row vector of non zero integers)
    GridSize;
    
    % Coordinates of the first vertex of the grid
    % (as a 1-by-2 row vector of double)
    GridOrigin;
    
    % Spacing between the vertices
    % (as a 1-by-2 row vector of double)
    GridSpacing;
end % end properties


%% Constructor
methods
    function obj = BSplineTransformModel2D(varargin)
        % Constructor for BSplineTransformModel2D class
        %
        % T = BSplineTransformModel2D();
        % Creates a new transform initialized with default values
        %
        % T = BSplineTransformModel2D(GRIDSIZE, GRIDSPACING, GRIDORIGIN);
        % Creates a new transform by specifying the grid parameters.
        %
        
        if nargin == 0
            % Initialization with default values
            nd = 2;
            obj.GridSize       = ones(1, nd);
            obj.GridSpacing    = ones(1, nd);
            obj.GridOrigin     = zeros(1, nd);
            initializeParameters();
                
        elseif nargin == 1
            % first argument is number of dimension
            var = varargin{1};
            if isscalar(var)
                nd = var;
                obj.GridSize       = ones(1, nd);
                obj.GridSpacing    = ones(1, nd);
                obj.GridOrigin     = zeros(1, nd);
                initializeParameters();
            end
            
        elseif nargin == 3
            obj.GridSize       = varargin{1};
            obj.GridSpacing    = varargin{2};
            obj.GridOrigin     = varargin{3};
            initializeParameters();
        end

        function initializeParameters()
            dim = obj.GridSize();
            np  = prod(dim) * length(dim);
            obj.Params = zeros(1, np);

            % initialize parameter names
            obj.ParamNames = cell(1, np);
            ind = 1;
            for iy = 1:obj.GridSize(2)
                for ix = 1:obj.GridSize(1)
                    obj.ParamNames{ind} = sprintf('vx_%d_%d', iy, ix);
                    ind = ind + 1;
                    obj.ParamNames{ind} = sprintf('vy_%d_%d', iy, ix);
                    ind = ind + 1;
                end
            end
        end

    end

end % end constructors


%% Methods specific to class
methods
    function res = subdivide(obj)
        %SUBDIVIDE Subdivide the transform grid.
        %
        % T2 = subdivide(T);
        
        % get 2D arrays of displacements in each direction
        dims = obj.GridSize;
        ux = reshape(obj.Params(1:2:end), dims);
        uy = reshape(obj.Params(2:2:end), dims);

        % subdivide first dimension
        dims2 = [dims(1)*2-1 dims(2)];
        ux2 = zeros(dims2);
        uy2 = zeros(dims2);
        ux2(1:2:end,:) = ux;
        uy2(1:2:end,:) = uy;
        ux2(2:2:end,:) = (ux(1:end-1,:) + ux(2:end,:)) / 2;
        uy2(2:2:end,:) = (uy(1:end-1,:) + uy(2:end,:)) / 2;
        ux = ux2;
        uy = uy2;
        dims = dims2;
        
        % subdivide second dimension
        dims2 = [dims(1) dims(2)*2-1];
        ux2 = zeros(dims2);
        uy2 = zeros(dims2);
        ux2(:,1:2:end) = ux;
        uy2(:,1:2:end) = uy;
        ux2(:,2:2:end) = (ux(:,1:end-1) + ux(:,2:end)) / 2;
        uy2(:,2:2:end) = (uy(:,1:end-1) + uy(:,2:end)) / 2;
        
        params = [ux2(:) uy2(:)]' / 2;
        params = params(:)';
        
        res = BSplineTransformModel2D(dims2, obj.GridSpacing/2, obj.GridOrigin);
        res.Params = params;
    end
    
    function drawVertexShifts(obj, varargin)
        % Draw the displacement associated to each vertex of the grid
        %
        % Example
        %    drawVertexShifts(T, 'g');
        %
        % See also
        %    drawGrid
        
        % get vertex array
        v = getGridVertices(obj);
        % get array of shifts
        shifts = getVertexShifts(obj);
        
        drawVector(v, shifts, varargin{:});
    end
    
    function drawGrid(obj)
        % Draw the grid used to defined the deformation
        % (Do not deform the grid)
        %
        % Example
        %    drawGrid(T);
        %
        % See also
        %    drawVertexShifts

        % create vertex array
        v = getGridVertices(obj);
        
        nv = size(v, 1);
        inds = reshape(1:nv, obj.GridSize);
        
        % edges in direction x
        ne1 = (obj.GridSize(2) - 1) * obj.GridSize(1);
        e1 = [reshape(inds(:, 1:end-1), [ne1 1]) reshape(inds(:, 2:end), [ne1 1])];
        
        % edges in direction y
        ne2 = obj.GridSize(2) * (obj.GridSize(1) - 1);
        e2 = [reshape(inds(1:end-1, :), [ne2 1]) reshape(inds(2:end, :), [ne2 1])];
        
        % create edge array
        e = cat(1, e1, e2);

        drawGraph(v, e);
    end
    
    function vertices = getGridVertices(obj)
        % Returns coordinates of grid vertices
        
        % base coordinates of grid vertices
        lx = (0:obj.GridSize(1) - 1) * obj.GridSpacing(1) + obj.GridOrigin(1);
        ly = (0:obj.GridSize(2) - 1) * obj.GridSpacing(2) + obj.GridOrigin(2);
        
        % create base mesh
        % (use reverse order to make vertices iterate in x order first)
        [x, y] = meshgrid(lx, ly);
        x = x';
        y = y';
        
        % create vertex array
        vertices = [x(:) y(:)];
    end
    
    function shifts = getVertexShifts(obj)
        % Returns shifts associated to each vertex as a N-by-2 array
        dx = reshape(obj.Params(1:2:end), obj.GridSize);
        dy = reshape(obj.Params(2:2:end), obj.GridSize);
        shifts = [dx(:) dy(:)];
    end
end


%% Modify or access the grid parameters
% the ix and iy parameters are the indices of the transform grid.
methods
    function ux = getUx(obj, ix, iy)
        ind = sub2ind(obj.GridSize, ix, iy) * 2 - 1;
        ux = obj.Params(ind);
    end
    
    function setUx(obj, ix, iy, ux)
        ind = sub2ind(obj.GridSize, ix, iy) * 2 - 1;
        obj.Params(ind) = ux;
    end
    
    function uy = getUy(obj, ix, iy)
        ind = sub2ind(obj.GridSize, ix, iy) * 2;
        uy = obj.Params(ind);
    end
    
    function setUy(obj, ix, iy, uy)
        ind = sub2ind(obj.GridSize, ix, iy) * 2;
        obj.Params(ind) = uy;
    end
end % end methods


%% Methods implementing the ParametricTransform interface
methods
    function jac = parametricJacobian(obj, x, varargin)
        % Compute parametric jacobian for a specific position.
        % 
        % jac = parametricJacobian(obj, x)
        % 
        % The result is a ND-by-NP array, where ND is the number of
        % dimension, and NP is the number of parameters.
        %
        % If x is a N-by-2 array, return result as a ND-by-NP-by-N array.
        %

        % extract coordinate of input point
        if isempty(varargin)
            y = x(:,2);
            x = x(:,1);
        else
            y = varargin{1};
        end

        % allocate result
        np = length(obj.Params);
        jac = zeros(2, np, length(x));
        dim = size(jac);
                
        % compute position wrt to the grid vertices (1-indexed)
        xg = (x - obj.GridOrigin(1)) / obj.GridSpacing(1) + 1;
        yg = (y - obj.GridOrigin(2)) / obj.GridSpacing(2) + 1;
        
        % coordinates within the unit tile
        xu = xg - floor(xg);
        yu = yg - floor(yg);
       
        baseFuns = {...
            @BSplines.beta3_0, ...
            @BSplines.beta3_1, ...
            @BSplines.beta3_2, ...
            @BSplines.beta3_3};
        
        % iteration on neighbor tiles 
        eval_i = zeros(size(xu));
        for i = -1:2
            % coordinates of neighbor grid vertex
            xv = floor(xg) + i;
            indOkX = xv >= 1 & xv <= obj.GridSize(1);

            % evaluate weight associated to grid vertex
            fun_i = baseFuns{i+2};
            eval_i(indOkX) = fun_i(xu(indOkX));
            
            for j = -1:2
                yv = floor(yg) + j;
                indOkY = yv >= 1 & yv <= obj.GridSize(2);

                % indices of points whose grid vertex is defined
                inds = indOkX & indOkY;
                
                % linear index of translation components
                indX = sub2ind(obj.GridSize, xv(inds), yv(inds)) * 2 - 1;
                
                % spline basis for y vertex
                fun_j = baseFuns{j+2};
                
                % evaluate weight associated to current grid vertex
                b = eval_i(inds) .* fun_j(yu(inds));
                
                % index of parameters
                indP = ones(size(indX));
                
                % update jacobian for grid vectors located around current
                % points
                jac(sub2ind(dim, indP, indX, find(inds))) = b;
                jac(sub2ind(dim, indP+1, indX+1, find(inds))) = b;
            end
        end
    end
end

%% Methods implementing the Transform interface
methods
    function point2 = transformPoint(obj, point)
        % Compute coordinates of transformed point
        
        % initialize coordinate of result
        point2 = point;
        
        % compute position wrt to the grid vertices (1-indexed)
        xg = (point(:, 1) - obj.GridOrigin(1)) / obj.GridSpacing(1) + 1;
        yg = (point(:, 2) - obj.GridOrigin(2)) / obj.GridSpacing(2) + 1;
        
        % coordinates within the unit tile
        xu = xg - floor(xg);
        yu = yg - floor(yg);
       
        baseFuns = {...
            @BSplines.beta3_0, ...
            @BSplines.beta3_1, ...
            @BSplines.beta3_2, ...
            @BSplines.beta3_3};
        
        % iteration on neighbor tiles 
        eval_i = zeros(size(xu));
        for i = -1:2
            % coordinates of neighbor grid vertex
            xv = floor(xg) + i;
            indOkX = xv >= 1 & xv <= obj.GridSize(1);

            % evaluate weight associated to grid vertex
            fun_i = baseFuns{i+2};
            eval_i(indOkX) = fun_i(xu(indOkX));
            
            for j = -1:2
                yv = floor(yg) + j;
                indOkY = yv >= 1 & yv <= obj.GridSize(2);

                % indices of points whose grid vertex is defined
                inds = indOkX & indOkY;
                
                % linear index of translation components
                indX = sub2ind(obj.GridSize, xv(inds), yv(inds)) * 2 - 1;
                
                % spline basis for y vertex
                fun_j = baseFuns{j+2};
                
                % evaluate weight associated to current grid vertex
                b = eval_i(inds) .* fun_j(yu(inds));
                
                % update coordinates of transformed points
                point2(inds,1) = point2(inds,1) + b .* obj.Params(indX)';
                point2(inds,2) = point2(inds,2) + b .* obj.Params(indX+1)';
            end
        end
    end
    
    function jac = jacobianMatrix(obj, point)
        % Jacobian matrix of the given point
        %
        %   JAC = getJacobian(TRANS, PT)
        %   where PT is a N-by-2 array of points, returns the spatial
        %   jacobian matrix of each point in the form of a 2-by-2-by-N
        %   array.
        %
        
        %% Constants
        
        % bspline basis functions and derivative functions
        baseFuns = {...
            @BSplines.beta3_0, ...
            @BSplines.beta3_1, ...
            @BSplines.beta3_2, ...
            @BSplines.beta3_3};
        
        derivFuns = {...
            @BSplines.beta3_0d, ...
            @BSplines.beta3_1d, ...
            @BSplines.beta3_2d, ...
            @BSplines.beta3_3d};

        
        %% Initializations
       
        % extract grid spacing for normalization
        deltaX = obj.GridSpacing(1);
        deltaY = obj.GridSpacing(2);
        
        % compute position of points wrt to grid vertices
        xg = (point(:, 1) - obj.GridOrigin(1)) / deltaX + 1;
        yg = (point(:, 2) - obj.GridOrigin(2)) / deltaY + 1;
        
        % initialize zeros translation vector
        nPts = length(xg);

        % coordinates within the unit tile
        xu = reshape(xg - floor(xg), [1 1 nPts]);
        yu = reshape(yg - floor(yg), [1 1 nPts]);       
        
        % allocate memory for storing result, and initialize to identity
        % matrix
        jac = zeros(2, 2, size(point, 1));
        jac(1, 1, :) = 1;
        jac(2, 2, :) = 1;
        
        % pre-allocate weights for vertex grids
        bx  = zeros(size(xu));
        bxd = zeros(size(xu));
        
        %% Iteration on neighbor tiles 
        for i = -1:2
            % x-coordinate of neighbor vertex
            xv  = floor(xg) + i;
            indOkX = xv >= 1 & xv <= obj.GridSize(1);
            
            % compute x-coefficients of bezier function and derivative
            bx(indOkX)  = baseFuns{i+2}(xu(indOkX));
            bxd(indOkX) = derivFuns{i+2}(xu(indOkX));
            
            for j = -1:2
                % y-coordinate of neighbor vertex
                yv = floor(yg) + j;
                indOkY = yv >= 1 & yv <= obj.GridSize(2);
                
                % indices of points whose grid vertex is defined
                inds = indOkX & indOkY;
                if all(~inds)
                    continue;
                end

                % linear index of translation components
                indX = sub2ind(obj.GridSize, xv(inds), yv(inds)) * 2 - 1;
                
                % translation vector of the current vertex
                dxv = reshape(obj.Params(indX), [1 1 length(indX)]);
                dyv = reshape(obj.Params(indX+1), [1 1 length(indX)]);
                
                % compute y-coefficients of bezier function and derivative
                by  = baseFuns{j+2}(yu(inds));
                byd = derivFuns{j+2}(yu(inds));

                % update jacobian matrix elements
                jac(1, 1, inds) = jac(1, 1, inds) + bxd(inds) .* by  .* dxv / deltaX;
                jac(1, 2, inds) = jac(1, 2, inds) + bx(inds)  .* byd .* dxv / deltaY;
                jac(2, 1, inds) = jac(2, 1, inds) + bxd(inds) .* by  .* dyv / deltaX;
                jac(2, 2, inds) = jac(2, 2, inds) + bx(inds)  .* byd .* dyv / deltaY;
            end
        end
    end

    function deriv = secondDerivatives(obj, point, indI, indJ)
        % Second derivatives for the given point(s)
        %
        % D2 = secondDerivatives(T, POINT, INDI, INDJ)
        % Return a M-by-2 array, with as many rows as the number of points.
        % First columns is the second derivative of the x-transform part,
        % and second column is the second derivative of the y-transform
        % part.
        
        %% Constants
        
        % bspline basis functions and derivative functions
        baseFuns = {...
            @BSplines.beta3_0, ...
            @BSplines.beta3_1, ...
            @BSplines.beta3_2, ...
            @BSplines.beta3_3};
        
        derivFuns = {...
            @BSplines.beta3_0d, ...
            @BSplines.beta3_1d, ...
            @BSplines.beta3_2d, ...
            @BSplines.beta3_3d};
        
        deriv2Funs = {...
            @BSplines.beta3_0s, ...
            @BSplines.beta3_1s, ...
            @BSplines.beta3_2s, ...
            @BSplines.beta3_3s};

        
        %% Initializations
       
        % extract grid spacing for normalization
        deltaX = obj.GridSpacing(1);
        deltaY = obj.GridSpacing(2);
        
        % compute position of points wrt to grid vertices
        xg = (point(:, 1) - obj.GridOrigin(1)) / deltaX + 1;
        yg = (point(:, 2) - obj.GridOrigin(2)) / deltaY + 1;
        
        % initialize zeros translation vector
        nPts = length(xg);

        % coordinates within the unit tile
        xu = reshape(xg - floor(xg), [nPts 1]);
        yu = reshape(yg - floor(yg), [nPts 1]);       
        
        % allocate memory for storing result
        deriv = zeros(size(point, 1), 2);
        
        % pre-allocate weights for vertex grids
        bx  = zeros(size(xu));
        bxd = zeros(size(xu));
        bxs = zeros(size(xu));
        
        %% Iteration on neighbor tiles 
        for i = -1:2
            % x-coordinate of neighbor vertex
            xv  = floor(xg) + i;
            indOkX = xv >= 1 & xv <= obj.GridSize(1);
            
            % compute x-coefficients of bezier function and derivative
            bx(indOkX)  = baseFuns{i+2}(xu(indOkX));
            bxd(indOkX) = derivFuns{i+2}(xu(indOkX));
            bxs(indOkX) = deriv2Funs{i+2}(xu(indOkX));
            
            for j = -1:2
                % y-coordinate of neighbor vertex
                yv = floor(yg) + j;
                indOkY = yv >= 1 & yv <= obj.GridSize(2);
                
                % indices of points whose grid vertex is defined
                inds = indOkX & indOkY;
                if all(~inds)
                    continue;
                end

                % linear index of translation components
                indX = sub2ind([obj.GridSize], xv(inds), yv(inds)) * 2 - 1;
                
                % translation vector of the current vertex
                dxv = reshape(obj.Params(indX), [length(indX) 1]);
                dyv = reshape(obj.Params(indX+1), [length(indX) 1]);
                
                % compute y-coefficients of spline function and derivative
                by  = baseFuns{j+2}(yu(inds));
                byd = derivFuns{j+2}(yu(inds));
                bys = deriv2Funs{j+2}(yu(inds));

                % update second derivatives elements
                if indI == 1 && indJ == 1
                    deriv(inds,1) = deriv(inds,1) + (bxs(inds) .* by  .* dxv) / (deltaX^2);
                    deriv(inds,2) = deriv(inds,2) + (bxs(inds) .* by  .* dyv) / (deltaX^2);
                    
                elseif indI == 2 && indJ == 2
                    deriv(inds,1) = deriv(inds,1) + (bx(inds)  .* bys .* dxv) / (deltaY^2);
                    deriv(inds,2) = deriv(inds,2) + (bx(inds)  .* bys .* dyv) / (deltaY^2);
                    
                elseif (indI == 1 && indJ == 2) || (indI == 2 && indJ == 1)
                    deriv(inds,1) = deriv(inds,1) + (bxd(inds) .* byd .* dxv) / (deltaX*deltaY);
                    deriv(inds,2) = deriv(inds,2) + (bxd(inds) .* byd .* dyv) / (deltaX*deltaY);

                else
                    error('indI and indJ should be between 1 and 2');
                end
            end
        end
        
    end % secondDerivatives

    function lap = curvatureOperator(obj, point)
        % Compute curvature operator at given position(s)
        %
        %   LAP = getLaplacian(TRANS, PT)
        %   where PT is a N-by-2 array of points, returns the laplacian of
        %   each point in the form of a 2-by-2-by-N array.
        %
        
        % compute second derivatives (each array is Npts-by-2
        dx2 = secondDerivatives(obj, point, 1, 1);
        dy2 = secondDerivatives(obj, point, 2, 2);
        
        % compute curvature operator
        lap = sum(dx2, 2).^2 + sum(dy2, 2).^2;
    end

    function dim = getDimension(obj) %#ok<MANU>
        dim = 2;
    end
end


%% Serialization methods
methods
    function str = toStruct(obj)
        % Converts to a structure to facilitate serialization
        str = struct('Type', 'BSplineTransformModel2D', ...
            'GridSize', obj.GridSize, ...
            'GridSpacing', obj.GridSpacing, ...
            'GridOrigin', obj.GridOrigin, ...
            'Parameters', obj.Params);
    end
end
methods (Static)
    function transfo = fromStruct(str)
        % Creates a new instance from a structure
        transfo = BSplineTransformModel2D(str.GridSize, str.GridSpacing, str.GridOrigin);
        transfo.Params = str.Parameters;
    end
end

end % end classdef

