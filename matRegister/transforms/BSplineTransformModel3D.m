classdef BSplineTransformModel3D < ParametricTransform
%BSPLINETRANSFORMMODEL3D Cubic B-Spline Transform model in 3D.
%
%   Class BSplineTransformModel3D
%
%   Grid is composed of M-by-N-by-P vertices, with M number of rows, N number
%   of columns and P number of planes. Iteration along x direction first.
%   Parameters correspond to shift vector associated to each vertex:
%   [vx111 vy111 vx211 vy211 ... vxIJK vyIJK ... vxMNP vyMNP]
%
%   Example
%   BSplineTransformModel3D
%
%   See also
%     BSplineTransformModel2D
%

% ------
% Author: David Legland
% e-mail: david.legland@inra.fr
% Created: 2018-08-09,    using Matlab 9.4.0.813654 (R2018a)
% Copyright 2018 INRA - BIA-BIBS.


%% Properties
properties
    % Number of vertices of the grid in each direction, in XYZ order.
    % (as a 1-by-3 row vector of non zero integers)
    GridSize;
    
    % Coordinates of the first vertex of the grid
    % (as a 1-by-3 row vector of double)
    GridOrigin;
    
    % Spacing between the vertices
    % (as a 1-by-3 row vector of double)
    GridSpacing;
end % end properties


%% Constructor
methods
    function obj = BSplineTransformModel3D(varargin)
        % Constructor for BSplineTransformModel3D class
        %
        % T = BSplineTransformModel3D();
        % Creates a new transform initialized with default values
        %
        % T = BSplineTransformModel3D(GRIDSIZE, GRIDSPACING, GRIDORIGIN);
        % Creates a new transform by specifying the grid parameters.
        %
        
        if nargin == 0
            % Initialization with default values
            nd = 3;
            obj.GridSize       = ones(1, nd);
            obj.GridSpacing    = ones(1, nd);
            obj.GridOrigin     = zeros(1, nd);
            initializeParameters();
                
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
            for iz = 1:obj.GridSize(3)
                for iy = 1:obj.GridSize(2)
                    for ix = 1:obj.GridSize(1)
                        obj.ParamNames{ind} = sprintf('vx_%d_%d_%d', ix, iy, iz);
                        ind = ind + 1;
                        obj.ParamNames{ind} = sprintf('vy_%d_%d_%d', ix, iy, iz);
                        ind = ind + 1;
                        obj.ParamNames{ind} = sprintf('vz_%d_%d_%d', ix, iy, iz);
                        ind = ind + 1;
                    end
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
        ux = reshape(obj.Params(1:3:end), dims);
        uy = reshape(obj.Params(2:3:end), dims);
        uz = reshape(obj.Params(3:3:end), dims);

        % subdivide first dimension
        dims2 = [dims(1)*2-1 dims(2) dims(3)];
        ux2 = zeros(dims2);
        uy2 = zeros(dims2);
        uz2 = zeros(dims2);
        ux2(1:2:end,:,:) = ux;
        uy2(1:2:end,:,:) = uy;
        uz2(1:2:end,:,:) = uz;
        ux2(2:2:end,:,:) = (ux(1:end-1,:,:) + ux(2:end,:,:)) / 2;
        uy2(2:2:end,:,:) = (uy(1:end-1,:,:) + uy(2:end,:,:)) / 2;
        uz2(2:2:end,:,:) = (uz(1:end-1,:,:) + uz(2:end,:,:)) / 2;
        ux = ux2;
        uy = uy2;
        uz = uz2;
        dims = dims2;
        
        % subdivide second dimension
        dims2 = [dims(1) dims(2)*2-1 dims(3)];
        ux2 = zeros(dims2);
        uy2 = zeros(dims2);
        uz2 = zeros(dims2);
        ux2(:,1:2:end,:) = ux;
        uy2(:,1:2:end,:) = uy;
        uz2(:,1:2:end,:) = uz;
        ux2(:,2:2:end,:) = (ux(:,1:end-1,:) + ux(:,2:end,:)) / 2;
        uy2(:,2:2:end,:) = (uy(:,1:end-1,:) + uy(:,2:end,:)) / 2;
        uz2(:,2:2:end,:) = (uz(:,1:end-1,:) + uz(:,2:end,:)) / 2;
        ux = ux2;
        uy = uy2;
        uz = uz2;
        dims = dims2;
        
        % subdivide third dimension
        dims2 = [dims(1) dims(2) dims(3)*2-1];
        ux2 = zeros(dims2);
        uy2 = zeros(dims2);
        uz2 = zeros(dims2);
        ux2(:,:,1:2:end) = ux;
        uy2(:,:,1:2:end) = uy;
        uz2(:,:,1:2:end) = uz;
        ux2(:,:,2:2:end) = (ux(:,:,1:end-1) + ux(:,:,2:end)) / 2;
        uy2(:,:,2:2:end) = (uy(:,:,1:end-1) + uy(:,:,2:end)) / 2;
        uz2(:,:,2:2:end) = (uz(:,:,1:end-1) + uz(:,:,2:end)) / 2;

        params = [ux2(:) uy2(:) uz2(:)]' / 2;
        params = params(:)';
        
        res = BSplineTransformModel3D(dims2, obj.GridSpacing/2, obj.GridOrigin);
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
        
        drawVector3d(v, shifts, varargin{:});
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
        
        nv = prod(obj.GridSize);
        inds = reshape(1:nv, obj.GridSize);
        
        nX = obj.GridSize(1);
        nY = obj.GridSize(2);
        nZ = obj.GridSize(3);
        
        % edges in direction x
        ne1 = (nX - 1) * nY * nZ;
        e1 = [reshape(inds(1:end-1, :, :), [ne1 1]) reshape(inds(2:end, :, :), [ne1 1])];
        
        % edges in direction y
        ne2 = nX * (nY - 1) * nZ;
        e2 = [reshape(inds(:, 1:end-1, :), [ne2 1]) reshape(inds(:, 2:end, :), [ne2 1])];
        
        % edges in direction z
        ne3 = nX * nY * (nZ - 1);
        e3 = [reshape(inds(:, :, 1:end-1), [ne3 1]) reshape(inds(:, :, 2:end), [ne3 1])];
        
        % create edge array
        e = cat(1, e1, e2, e3);

        drawGraph(v, e);
    end
    
    function vertices = getGridVertices(obj)
        % Returns coordinates of grid vertices
        
        % base coordinates of grid vertices
        lx = (0:obj.GridSize(1) - 1) * obj.GridSpacing(1) + obj.GridOrigin(1);
        ly = (0:obj.GridSize(2) - 1) * obj.GridSpacing(2) + obj.GridOrigin(2);
        lz = (0:obj.GridSize(3) - 1) * obj.GridSpacing(3) + obj.GridOrigin(3);
        
        % create base mesh
        % (use reverse order to make vertices iterate in x order first)
        [x, y, z] = meshgrid(lx, ly, lz);
        x = permute(x, [2 1 3]);
        y = permute(y, [2 1 3]);
        z = permute(z, [2 1 3]);
        
        % create vertex array
        vertices = [x(:) y(:) z(:)];
    end
    
    function shifts = getVertexShifts(obj)
        % Returns shifts associated to each vertex as a N-by-3 array
        dx = reshape(obj.Params(1:3:end), obj.GridSize);
        dy = reshape(obj.Params(2:3:end), obj.GridSize);
        dz = reshape(obj.Params(3:3:end), obj.GridSize);
        shifts = [dx(:) dy(:) dz(:)];
    end
end


%% Modify or access the grid parameters
% the ix and iy parameters are the indices of the transform grid.
methods
    function ux = getUx(obj, ix, iy)
        ind = sub2ind(obj.GridSize, ix, iy, iz) * 3 - 2;
        ux = obj.Params(ind);
    end
    
    function setUx(obj, ix, iy, iz, ux)
        ind = sub2ind(obj.GridSize, ix, iy, iz) * 3 - 2;
        obj.Params(ind) = ux;
    end
    
    function uy = getUy(obj, ix, iy, iz)
        ind = sub2ind(obj.GridSize, ix, iy, iz) * 3 - 1;
        uy = obj.Params(ind);
    end
    
    function setUy(obj, ix, iy, iz, uy)
        ind = sub2ind(obj.GridSize, ix, iy, iz) * 3 - 1;
        obj.Params(ind) = uy;
    end
    
    function uz = getUz(obj, ix, iy, iz)
        ind = sub2ind(obj.GridSize, ix, iy, iz) * 3;
        uz = obj.Params(ind);
    end
    
    function setUz(obj, ix, iy, iz, uz)
        ind = sub2ind(obj.GridSize, ix, iy, iz) * 3;
        obj.Params(ind) = uz;
    end
end % end methods


%% Methods implementing the ParametricTransform interface
methods
    function jac = parametricJacobian(obj, x, varargin)
        % Computes parametric jacobian for a specific position
        % 
        % jac = getParametricJacobian(obj, x)
        % 
        % The result is a ND-by-NP array, where ND is the number of
        % dimension, and NP is the number of parameters.
        %
        % If x is a N-by-3 array, return result as a ND-by-NP-by-N array.
        %

        % extract coordinate of input point
        if isempty(varargin)
            y = x(:,2);
            z = x(:,3);
            x = x(:,1);
        else
            y = varargin{1};
            z = varargin{2};
        end

        % allocate result
        np = length(obj.Params);
        jac = zeros(3, np, length(x));
        dim = size(jac);
                
        % compute position wrt to the grid vertices (1-indexed)
        xg = (x - obj.GridOrigin(1)) / obj.GridSpacing(1) + 1;
        yg = (y - obj.GridOrigin(2)) / obj.GridSpacing(2) + 1;
        zg = (z - obj.GridOrigin(3)) / obj.GridSpacing(3) + 1;
        
        % coordinates within the unit tile
        xu = xg - floor(xg);
        yu = yg - floor(yg);
        zu = zg - floor(zg);
       
        baseFuns = {...
            @BSplines.beta3_0, ...
            @BSplines.beta3_1, ...
            @BSplines.beta3_2, ...
            @BSplines.beta3_3};
        
        
        % TODO: test this function
        % iteration on neighbor tiles
        eval_j = zeros(size(xu));
        eval_k = zeros(size(xu));
        for k = -1:2
            % coordinates of neighbor grid vertex
            zv = floor(zg) + k;
            indOkZ = zv >= 1 & zv <= obj.GridSize(3);

            % evaluate weight associated to grid vertex
            fun_k = baseFuns{k+2};
            eval_k(indOkZ) = fun_k(zu(indOkZ));
            
            for j = -1:2
                % coordinates of neighbor grid vertex
                yv = floor(yg) + j;
                indOkY = yv >= 1 & yv <= obj.GridSize(2);
            
                % evaluate weight associated to grid vertex
                fun_j = baseFuns{j+2};
                eval_j(indOkY) = fun_j(yu(indOkY));
                
                for i = -1:2
%                     fprintf('%d,%d,%d\n', i, j, k);
                    
                    % coordinates of neighbor grid vertex
                    xv = floor(xg) + i;
                    indOkX = xv >= 1 & xv <= obj.GridSize(1);
                
                    % indices of points whose grid vertex is defined
                    inds = indOkX & indOkY & indOkZ;
                    
                    % linear index of translation components
                    indX = sub2ind(obj.GridSize, xv(inds), yv(inds), zv(inds)) * 3 - 2;
                    
                    % spline basis for x vertex
                    fun_i = baseFuns{i+2};
                    
                    % evaluate weight associated to current grid vertex
                    b = fun_i(xu(inds)) .* eval_j(inds) .* eval_k(inds);
                    
                    % update jacobian for grid vectors located around current
                    % points
                    jac(sub2ind(dim, indP, indX, find(inds))) = b;
                    jac(sub2ind(dim, indP+1, indX+1, find(inds))) = b;
                    jac(sub2ind(dim, indP+2, indX+2, find(inds))) = b;
                end
            end
        end
    end
    
    function transfo = clone(obj)
        transfo = BSplineTransformModel3D(obj.GridSize, obj.GridSpacing, obj.GridOrigin);
        transfo.Params = obj.Params;
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
        zg = (point(:, 3) - obj.GridOrigin(3)) / obj.GridSpacing(3) + 1;
        
        % coordinates within the unit tile
        xu = xg - floor(xg);
        yu = yg - floor(yg);
        zu = zg - floor(zg);
       
        baseFuns = {...
            @BSplines.beta3_0, ...
            @BSplines.beta3_1, ...
            @BSplines.beta3_2, ...
            @BSplines.beta3_3};
        
        % iteration on neighbor tiles 
        eval_j = zeros(size(xu));
        eval_k = zeros(size(xu));
        for k = -1:2
            % coordinates of neighbor grid vertex
            zv = floor(zg) + k;
            indOkZ = zv >= 1 & zv <= obj.GridSize(3);

            % evaluate weight associated to grid vertex
            fun_k = baseFuns{k+2};
            eval_k(indOkZ) = fun_k(zu(indOkZ));
            
            for j = -1:2
                % coordinates of neighbor grid vertex
                yv = floor(yg) + j;
                indOkY = yv >= 1 & yv <= obj.GridSize(2);
            
                % evaluate weight associated to grid vertex
                fun_j = baseFuns{j+2};
                eval_j(indOkY) = fun_j(yu(indOkY));
                
                for i = -1:2
%                     fprintf('%d,%d,%d\n', i, j, k);
                    
                    % coordinates of neighbor grid vertex
                    xv = floor(xg) + i;
                    indOkX = xv >= 1 & xv <= obj.GridSize(1);
                
                    % indices of points whose grid vertex is defined
                    inds = indOkX & indOkY & indOkZ;
                    
                    % linear index of translation components
                    indX = sub2ind(obj.GridSize, xv(inds), yv(inds), zv(inds)) * 3 - 2;
                    
                    % spline basis for x vertex
                    fun_i = baseFuns{i+2};
                    
                    % evaluate weight associated to current grid vertex
                    b = fun_i(xu(inds)) .* eval_j(inds) .* eval_k(inds);
                    
                    % update coordinates of transformed points
                    point2(inds,1) = point2(inds,1) + b .* obj.Params(indX)';
                    point2(inds,2) = point2(inds,2) + b .* obj.Params(indX+1)';
                    point2(inds,3) = point2(inds,3) + b .* obj.Params(indX+2)';
                end
            end
        end
    end
    
    function jac = jacobianMatrix(obj, point)
        % Compute Jacobian matrix at the given point.
        %
        %   JAC = jacobianMatrix(TRANSFO, PT)
        %   where PT is a N-by-3 array of points, returns the spatial
        %   jacobian matrix of each point in the form of a 3-by-3-by-N
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
        deltaZ = obj.GridSpacing(3);
        
        % compute position of points wrt to grid vertices
        xg = (point(:, 1) - obj.GridOrigin(1)) / deltaX + 1;
        yg = (point(:, 2) - obj.GridOrigin(2)) / deltaY + 1;
        zg = (point(:, 3) - obj.GridOrigin(3)) / deltaZ + 1;
        
        % initialize zeros translation vector
        nPts = length(xg);

        % coordinates within the unit tile
        xu = reshape(xg - floor(xg), [1 1 nPts]);
        yu = reshape(yg - floor(yg), [1 1 nPts]);       
        zu = reshape(zg - floor(zg), [1 1 nPts]);       
        
        % allocate memory for storing result, and initialize to identity
        % matrix
        jac = zeros(3, 3, size(point, 1));
        jac(1, 1, :) = 1;
        jac(2, 2, :) = 1;
        jac(3, 3, :) = 1;
        
        % pre-allocate weights for vertex grids
        bz  = zeros(size(zu));
        bzd = zeros(size(zu));
        by  = zeros(size(yu));
        byd = zeros(size(yu));
                
        %% Iteration on neighbor tiles
        for k = -1:2
            % y-coordinate of neighbor vertex
            zv = floor(zg) + k;
            indOkZ = zv >= 1 & zv <= obj.GridSize(3);
            
            % compute z-coefficients of bezier function and derivative
            bz(indOkZ)  = baseFuns{k+2}(zu(indOkZ));
            bzd(indOkZ) = derivFuns{k+2}(zu(indOkZ));
            
            for j = -1:2
                % y-coordinate of neighbor vertex
                yv = floor(yg) + j;
                indOkY = yv >= 1 & yv <= obj.GridSize(2);
                
                % compute y-coefficients of bezier function and derivative
                by(indOkY)  = baseFuns{j+2}(yu(indOkY));
                byd(indOkY) = derivFuns{j+2}(yu(indOkY));
                
                for i = -1:2
                    % x-coordinate of neighbor vertex
                    xv  = floor(xg) + i;
                    indOkX = xv >= 1 & xv <= obj.GridSize(1);
                    
                    % indices of points whose grid vertex is defined
                    inds = indOkX & indOkY & indOkZ;
                    if all(~inds)
                        continue;
                    end
                    
                    % linear index of translation components
                    indX = sub2ind(obj.GridSize, xv(inds), yv(inds), zv(inds)) * 3 - 2;
                    
                    % translation vector of the current vertex
                    dxv = reshape(obj.Params(indX),   [1 1 length(indX)]);
                    dyv = reshape(obj.Params(indX+1), [1 1 length(indX)]);
                    dzv = reshape(obj.Params(indX+2), [1 1 length(indX)]);
                    
                    % compute x-coefficients of bezier function and derivative
                    bx  = baseFuns{i+2}(xu(inds));
                    bxd = derivFuns{i+2}(xu(inds));
                    
                    % update elements of the 3-by-3 jacobian matrices
                    jac(1, 1, inds) = jac(1, 1, inds) + bxd .* by(inds) .* bz(inds) .* dxv / deltaX;
                    jac(1, 2, inds) = jac(1, 2, inds) + bx .* byd(inds) .* bz(inds) .* dxv / deltaY;
                    jac(1, 3, inds) = jac(1, 3, inds) + bx .* by(inds) .* bzd(inds) .* dxv / deltaZ;
                    jac(2, 1, inds) = jac(2, 1, inds) + bxd .* by(inds) .* bz(inds) .* dyv / deltaX;
                    jac(2, 2, inds) = jac(2, 2, inds) + bx .* byd(inds) .* bz(inds) .* dyv / deltaY;
                    jac(2, 3, inds) = jac(2, 3, inds) + bx .* by(inds) .* bzd(inds) .* dyv / deltaZ;
                    jac(3, 1, inds) = jac(3, 1, inds) + bxd .* by(inds) .* bz(inds) .* dzv / deltaX;
                    jac(3, 2, inds) = jac(3, 2, inds) + bx .* byd(inds) .* bz(inds) .* dzv / deltaY;
                    jac(3, 3, inds) = jac(3, 3, inds) + bx .* by(inds) .* bzd(inds) .* dzv / deltaZ;
                end
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
        deltaZ = obj.GridSpacing(3);
        
        % compute position of points wrt to grid vertices
        xg = (point(:, 1) - obj.GridOrigin(1)) / deltaX + 1;
        yg = (point(:, 2) - obj.GridOrigin(2)) / deltaY + 1;
        zg = (point(:, 3) - obj.GridOrigin(3)) / deltaZ + 1;
        
        % initialize zeros translation vector
        nPts = length(xg);

        % coordinates within the unit tile
        xu = reshape(xg - floor(xg), [nPts 1]);
        yu = reshape(yg - floor(yg), [nPts 1]);
        zu = reshape(zg - floor(zg), [nPts 1]);
        
        % allocate memory for storing result
        deriv = zeros(size(point, 1), 3);
        
        % pre-allocate weights for vertex grids
        bz  = zeros(size(xu));
        bzd = zeros(size(xu));
        bzs = zeros(size(xu));
        by  = zeros(size(xu));
        byd = zeros(size(xu));
        bys = zeros(size(xu));
        
        %% Iteration on neighbor tiles 
        for k = -1:2
            % z-coordinate of neighbor vertex
            zv = floor(zg) + k;
            indOkZ = zv >= 1 & zv <= obj.GridSize(3);
            
            % compute z-coefficients of bezier function and derivative
            bz(indOkZ)  = baseFuns{k+2}(zu(indOkZ));
            bzd(indOkZ) = derivFuns{k+2}(zu(indOkZ));
            bzs(indOkZ) = deriv2Funs{k+2}(zu(indOkZ));
            
            for j = -1:2
                % y-coordinate of neighbor vertex
                yv = floor(yg) + j;
                indOkY = yv >= 1 & yv <= obj.GridSize(2);
                
                % compute y-coefficients of bezier function and derivative
                by(indOkY)  = baseFuns{j+2}(yu(indOkY));
                byd(indOkY) = derivFuns{j+2}(yu(indOkY));
                bys(indOkY) = deriv2Funs{j+2}(yu(indOkY));
                
                for i = -1:2
                    % x-coordinate of neighbor vertex
                    xv  = floor(xg) + i;
                    indOkX = xv >= 1 & xv <= obj.GridSize(1);
                    
                    % indices of points whose grid vertex is defined
                    inds = indOkX & indOkY & indOkZ;
                    if all(~inds)
                        continue;
                    end
                    
                    % linear index of translation components
                    indX = sub2ind([obj.GridSize], xv(inds), yv(inds), zv(inds)) * 3 - 2;
                    
                    % translation vector of the current vertex
                    dxv = reshape(obj.Params(indX),   [length(indX) 1]);
                    dyv = reshape(obj.Params(indX+1), [length(indX) 1]);
                    dzv = reshape(obj.Params(indX+2), [length(indX) 1]);
                    
                    % compute x-coefficients of spline function and derivative
                    bx  = baseFuns{i+2}(xu(inds));
                    bxd = derivFuns{i+2}(xu(inds));
                    bxs = deriv2Funs{i+2}(xu(inds));
                    
                    % update second derivatives elements
                    if indI == 1 && indJ == 1
                        deriv(inds,1) = deriv(inds,1) + (bxs .* by(inds) .* bz(inds) .* dxv) / (deltaX^2);
                        deriv(inds,2) = deriv(inds,2) + (bxs .* by(inds) .* bz(inds) .* dyv) / (deltaX^2);
                        deriv(inds,3) = deriv(inds,3) + (bxs .* by(inds) .* bz(inds) .* dzv) / (deltaX^2);
                        
                    elseif indI == 2 && indJ == 2
                        deriv(inds,1) = deriv(inds,1) + (bx .* bys(inds) .* bz(inds) .* dxv) / (deltaY^2);
                        deriv(inds,2) = deriv(inds,2) + (bx .* bys(inds) .* bz(inds) .* dyv) / (deltaY^2);
                        deriv(inds,3) = deriv(inds,3) + (bx .* bys(inds) .* bz(inds) .* dzv) / (deltaY^2);
                        
                    elseif indI == 3 && indJ == 3
                        deriv(inds,1) = deriv(inds,1) + (bx .* by(inds) .* bzs(inds) .* dxv) / (deltaZ^2);
                        deriv(inds,2) = deriv(inds,2) + (bx .* by(inds) .* bzs(inds) .* dyv) / (deltaZ^2);
                        deriv(inds,3) = deriv(inds,3) + (bx .* by(inds) .* bzs(inds) .* dzv) / (deltaZ^2);
                        
                    elseif (indI == 1 && indJ == 2) || (indI == 2 && indJ == 1)
                        deriv(inds,1) = deriv(inds,1) + (bxd .* byd(inds) .* bz(inds) .* dxv) / (deltaX*deltaY);
                        deriv(inds,2) = deriv(inds,2) + (bxd .* byd(inds) .* bz(inds) .* dyv) / (deltaX*deltaY);
                        deriv(inds,3) = deriv(inds,3) + (bxd .* byd(inds) .* bz(inds) .* dzv) / (deltaX*deltaY);
                        
                    elseif (indI == 1 && indJ == 3) || (indI == 3 && indJ == 1)
                        deriv(inds,1) = deriv(inds,1) + (bxd .* by(inds) .* bzd(inds) .* dxv) / (deltaX*deltaZ);
                        deriv(inds,2) = deriv(inds,2) + (bxd .* by(inds) .* bzd(inds) .* dyv) / (deltaX*deltaZ);
                        deriv(inds,3) = deriv(inds,3) + (bxd .* by(inds) .* bzd(inds) .* dzv) / (deltaX*deltaZ);
                        
                    elseif (indI == 2 && indJ == 3) || (indI == 3 && indJ == 2)
                        deriv(inds,1) = deriv(inds,1) + (bx .* byd(inds) .* bzd(inds) .* dxv) / (deltaY*deltaZ);
                        deriv(inds,2) = deriv(inds,2) + (bx .* byd(inds) .* bzd(inds) .* dyv) / (deltaY*deltaZ);
                        deriv(inds,3) = deriv(inds,3) + (bx .* byd(inds) .* bzd(inds) .* dzv) / (deltaY*deltaZ);
                        
                    else
                        error('indI and indJ should be between 1 and 2');
                    end
                end
            end
        end
        
    end % secondDerivatives

    function lap = curvatureOperator(obj, point)
        % Compute curvature operator at given position(s)
        %
        %   LAP = curvatureOperator(TRANS, PT)
        %   where PT is a N-by-3 array of points, returns the laplacian of
        %   each point in the form of a 3-by-3-by-N array.
        %
        
        % compute second derivatives (each array is Npts-by-2
        dx2 = secondDerivatives(obj, point, 1, 1);
        dy2 = secondDerivatives(obj, point, 2, 2);
        dz2 = secondDerivatives(obj, point, 3, 3);
        
        % compute curvature operator
        lap = sum(dx2, 2).^2 + sum(dy2, 2).^2 + sum(dz2, 2).^2;
    end

    function dim = getDimension(obj) %#ok<MANU>
        dim = 3;
    end
end


%% Serialization methods
methods
    function str = toStruct(obj)
        % Converts to a structure to facilitate serialization
        str = struct('Type', 'BSplineTransformModel3D', ...
            'GridSize', obj.GridSize, ...
            'GridSpacing', obj.GridSpacing, ...
            'GridOrigin', obj.GridOrigin, ...
            'Parameters', obj.Params);
    end
end
methods (Static)
    function transfo = fromStruct(str)
        % Creates a new instance from a structure
        transfo = BSplineTransformModel3D(str.GridSize, str.GridSpacing, str.GridOrigin);
        transfo.Params = str.Parameters;
    end
end

end % end classdef


