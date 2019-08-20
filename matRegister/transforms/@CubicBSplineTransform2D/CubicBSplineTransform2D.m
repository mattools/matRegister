classdef CubicBSplineTransform2D < ParametricTransform
%CUBICBSPLINETRANSFORM2D Cubic Spline Transform model in 2D
%
%   output = CubicBSplineTransform2D(input)
%
%   Example
%   CubicBSplineTransform2D
%
%   See also
%

% ------
% Author: David Legland
% e-mail: david.legland@inra.fr
% Created: 2011-02-16,    using Matlab 7.9.0.529 (R2009b)
% Copyright 2011 INRA - Cepia Software Platform.

%% Properties
properties
    % Number of vertices of the grid in each direction
    GridSize;
    
    % Coordinates of the first vertex of the grid
    GridOrigin;
    
    % Spacing between the vertices
    GridSpacing;
end 

%% Constructor
methods
    function obj = CubicBSplineTransform2D(varargin)
        % Create a new CubicBSplineTransform2D
        %
        % TRANSFO = CubicBSplineTransform2D();
        % Creates a new transfo initialized with default values
        %
        % TRANSFO = CubicBSplineTransform2D(GRIDSIZE, GRIDSPACING, GRIDORIGIN);
        % Creates a new transfo by specifying the grid parameters.
                
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
            dim = obj.GridSize;
            np  = prod(dim) * length(dim);
            obj.Params = zeros(1, np);
        end
    end % constructor 

end % construction function

%% General methods
methods

    function dim = getDimension(obj) %#ok<MANU>
        dim = 2;
    end

    function [point2, isInside] = transformPoint(obj, point)
        % Compute corodinates of transformed point
        
        % compute position wrt to the grid vertices
        xg = (point(:, 1) - obj.GridOrigin(1)) / obj.GridSpacing(1) + 1;
        yg = (point(:, 2) - obj.GridOrigin(2)) / obj.GridSpacing(2) + 1;
        
        % compute indices of points located within interpolation area
        isInsideX   = xg >= 2 & xg < obj.GridSize(1)-1;
        isInsideY   = yg >= 2 & yg < obj.GridSize(2)-1;
        isInside    = isInsideX & isInsideY;
        
        % select valid points
        xg = xg(isInside);
        yg = yg(isInside);
        
        % compute indices in linear indexing
        dimGrid = obj.GridSize;
        dimXY = dimGrid(1) * dimGrid(2);

        % coordinates within the unit tile
        xu = xg - floor(xg);
        yu = yg - floor(yg);
       
        % initialize zeros translation vector
        dx = zeros(length(xg), 1);
        dy = zeros(length(xg), 1);
        
        baseFuns = {...
            @BSplines.beta3_0, ...
            @BSplines.beta3_1, ...
            @BSplines.beta3_2, ...
            @BSplines.beta3_3};
        
        % iteration on each tile of the grid
        for i = -1:2
            % coordinates of neighbor vertex
            xv = floor(xg) + i;
            
            fun_i = baseFuns{i+2};
            eval_i = fun_i(xu);
            
            for j = -1:2
                yv = floor(yg) + j;
                
                % linear index of translation components
                indX = sub2ind([obj.GridSize], xv, yv);
                indY = sub2ind([obj.GridSize], xv, yv) + dimXY;
                
                % translation vector of the current vertex
                dxv = obj.Params(indX)';
                dyv = obj.Params(indY)';
                
                fun_j = baseFuns{j+2};
                
                % update total translation component
                b  = eval_i .* fun_j(yu);
                dx = dx + b.*dxv; 
                dy = dy + b.*dyv; 
            end
        end
        
        % update coordinates of transformed point
        point2 = point;
        point2(isInside, 1) = point(isInside, 1) + dx;
        point2(isInside, 2) = point(isInside, 2) + dy;
                
    end
    
    function ux = getUx(obj, x, y)
        ind = sub2ind([obj.GridSize 2], x, y, 1);
        ux = obj.Params(ind);
    end
    
    function setUx(obj, x, y, ux)
        ind = sub2ind([obj.GridSize 2], x, y, 1);
        obj.Params(ind) = ux;
    end
    
    function uy = getUy(obj, x, y)
        ind = sub2ind([obj.GridSize 2], x, y, 2);
        uy = obj.Params(ind);
    end
    
    function setUy(obj, x, y, uy)
        ind = sub2ind([obj.GridSize 2], x, y, 2);
        obj.Params(ind) = uy;
    end
    
    function drawGrid(obj)

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
        % get coordinates of grid vertices
        
        % base coordinates of grid vertices
        lx = (0:obj.GridSize(1) - 1) * obj.GridSpacing(1) + obj.GridOrigin(1);
        ly = (0:obj.GridSize(2) - 1) * obj.GridSpacing(2) + obj.GridOrigin(2);
        
        % create base mesh
        [x, y] = meshgrid(lx, ly);
        
        % add grid shifts
        x = x' + reshape(obj.Params(1:end/2), obj.GridSize);
        y = y' + reshape(obj.Params(end/2+1:end), obj.GridSize);
        
        % create vertex array
        vertices = [x(:) y(:)];
    end
    
    function transformVector(obj, varargin)
        error('MatRegister:UnimplementedMethod', ...
            'Method "%s" is not implemented for class "%s"', ...
            'transformVector', mfilename);
    end
    
    function jac = jacobianMatrix(obj, point)
        % Jacobian matrix of the given point
        %
        %   JAC = getJacobian(TRANS, PT)
        %   where PT is a N-by-2 array of points, returns the spatial
        %   jacobian matrix of each poitn in the form of a 2-by-2-by-N
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
        
        % compute indices of values within interpolation area
        isInsideX = xg >= 2 & xg < obj.GridSize(1)-1;
        isInsideY = yg >= 2 & yg < obj.GridSize(2)-1;
        inds = isInsideX & isInsideY;
        
        % keep only valid positions
        xg = xg(inds);
        yg = yg(inds);
        
        % initialize zeros translation vector
        nValid = length(xg);

        % coordinates within the unit tile
        xu = reshape(xg - floor(xg), [1 1 nValid]);
        yu = reshape(yg - floor(yg), [1 1 nValid]);       

        % compute indices in linear indexing
        dimXY = prod(obj.GridSize);
        
        % allocate memory for storing result, and initialize to identity
        % matrix
        jac = zeros(2, 2, size(point, 1));
        jac(1, 1, :) = 1;
        jac(2, 2, :) = 1;
        
        
        %% Iteration on neighbor tiles 
        
        for i = -1:2
            % x-coordinate of neighbor vertex
            xv  = floor(xg) + i;
            
            % compute x-coefficients of bezier function and derivative
            bx  = baseFuns{i+2}(xu);
            bxd = derivFuns{i+2}(xu);
            
            for j = -1:2
                % y-coordinate of neighbor vertex
                yv = floor(yg) + j;
                
                % linear index of translation components
                indX = sub2ind([obj.GridSize], xv, yv);
                indY = sub2ind([obj.GridSize], xv, yv) + dimXY;
                
                % translation vector of the current vertex
                dxv = reshape(obj.Params(indX), [1 1 length(inds)]);
                dyv = reshape(obj.Params(indY), [1 1 length(inds)]);
                
                % compute y-coefficients of bezier function and derivative
                by  = baseFuns{j+2}(yu);
                byd = derivFuns{j+2}(yu);

                % update jacobian matrix elements
                jac(1, 1, inds) = jac(1, 1, inds) + bxd .* by  .* dxv / deltaX;
                jac(1, 2, inds) = jac(1, 2, inds) + bx  .* byd .* dxv / deltaY;
                jac(2, 1, inds) = jac(2, 1, inds) + bxd .* by  .* dyv / deltaX;
                jac(2, 2, inds) = jac(2, 2, inds) + bx  .* byd .* dyv / deltaY;
            end
        end

    end
    
    function jac = parametricJacobian(obj, x, varargin)
        % Compute parametric jacobian for a specific position
        % The result is a ND-by-NP array, where ND is the number of
        % dimension, and NP is the number of parameters.
        
        % extract coordinate of input point
        if isempty(varargin)
            y = x(:,2);
            x = x(:,1);
        else
            y = varargin{1};
        end
        
        % compute position wrt to the grid vertices
        deltaX = obj.GridSpacing(1);
        deltaY = obj.GridSpacing(2);
        xg = (x - obj.GridOrigin(1)) / deltaX + 1;
        yg = (y - obj.GridOrigin(2)) / deltaY + 1;
        
        % compute indices of values within interpolation area
        isInsideX = xg >= 2 & xg < obj.GridSize(1)-1;
        isInsideY = yg >= 2 & yg < obj.GridSize(2)-1;
        isInside = isInsideX & isInsideY;
        inds = find(isInside);
        
        % keep only valid positions
        xg = xg(isInside);
        yg = yg(isInside);
        
        % initialize zeros translation vector
        nValid = length(xg);

        % pre-allocate result array
        nd = length(obj.GridSize);
        np = length(obj.Params);
        jac = zeros(nd, np, length(x));

        % if point is outside, return zeros matrix
        if ~isInside
            return;
        end
        
        % coordinates within the unit tile
        xu = reshape(xg - floor(xg), [1 1 nValid]);
        yu = reshape(yg - floor(yg), [1 1 nValid]);       
        
        dimGrid = obj.GridSize;
        dimX    = dimGrid(1);
        dimXY   = dimX * dimGrid(2);
        
        baseFuns = {...
            @BSplines.beta3_0, ...
            @BSplines.beta3_1, ...
            @BSplines.beta3_2, ...
            @BSplines.beta3_3};
        
        
        % pre-compute values of b-splines functions
        evals_i = zeros(nValid, 4);
        evals_j = zeros(nValid, 4);
        for i = 1:4
            fun_i = baseFuns{i};
            evals_i(:,i) = fun_i(xu);
            evals_j(:,i) = fun_i(yu);
        end
        
        % iteration on each tile of the grid
        for i = -1:2
            xv = floor(xg) + i;
            
            for j = -1:2
                % coordinates of neighbor vertex
                yv = floor(yg) + j;
                
                % linear index of translation components
                indX = (yv - 1) * dimX + xv;
                indY = indX + dimXY;
                                
                % update total translation component
                b = evals_i(:,i+2) .* evals_j(:,j+2);
                
                % update jacobian matrix (of size nd * np * nPts)
                jacInds = (indX - 1 + (inds - 1) * np) * 2 + 1;
                jac(jacInds) = b;
                jacInds = (indY - 1 + (inds - 1) * np) * 2 + 1;
                jac(jacInds + 1) = b;
%                 % equivalent to:
%                 iValid = ones(nValid, 1);
%                 jac(sub2ind(size(jac), 1*iValid, indX, inds)) = b;
%                 jac(sub2ind(size(jac), 2*iValid, indY, inds)) = b;
            end
        end
        
    end
    
   
end % general methods


%% I/O Methods
methods
    function writeToFile(obj, file)
        % Write transform parameter to the given file handle
        % Assumes file handle is an instance of FileWriter.
        %
        % Example
        %   F = fopen('transfo.txt', 'wt');
        %   fprintf(F, '#--- Transform Parameters ---');
        %   writeToFile(TRANSFO, F);
        %   fclose(F);
        %
        
        closeFile = false;
        if ischar(file)
            file = fopen(file, 'wt');
            closeFile = true;
        end
        
        nDims = 2;
        
        fprintf(file, 'TransformType = %s\n', class(obj));
        fprintf(file, 'TransformDimension = %d\n', nDims);
        
        nParams = length(obj.Params);
        fprintf(file, 'TransformParameterNumber = %d \n', nParams);
        
        pattern = ['TransformParameters =', repmat(' %g', 1, nParams) '\n'];
        fprintf(file, pattern, obj.Params);
        
        % some transform specific settings
        pattern = ['%s =' repmat(' %g', 1, nDims) '\n'];
        fprintf(file, pattern, 'TransformGridSize',     obj.GridSize);
        fprintf(file, pattern, 'TransformGridOrigin',   obj.GridOrigin);
        fprintf(file, pattern, 'TransformGridSpacing',  obj.GridSpacing);
        
        % close file
        if closeFile
            fclose(file);
        end
    end
end

methods (Static)
    function transfo = readFromFile(fileName)
        % Read transform from the given file name.
        % Returns a new instance of CubicBSplineTransform2D.
        %
        % Example
        %   TRANSFO = CubicBSplineTransform2D.readFromFile('transfo.txt');
        
        map = readPropertyFile(fileName);
        transfo = CubicBSplineTransform2D.createFromPropertyMap(map);
    end
    
    function transfo = createFromPropertyMap(map)
        % Create a new transform from a set of properties
        
        grSize  = map('TransformGridSize');
        grSize  = cellfun(@str2double, regexp(grSize, '\s*', 'split'));
        grSpac  = map('TransformGridSpacing');
        grSpac  = cellfun(@str2double, regexp(grSpac, '\s*', 'split'));
        grOrig  = map('TransformGridOrigin');
        grOrig  = cellfun(@str2double, regexp(grOrig, '\s*', 'split'));
        
        transfo = CubicBSplineTransform2D(grSize, grSpac, grOrig);
        
        
        nbParams = str2double(map('TransformParameterNumber'));
        
        trParams = map('TransformParameters');
        trParams= cellfun(@str2double, regexp(trParams, '\s*', 'split'));
        
        if nbParams ~= length(trParams)
            error('Wrong number of parameters');
        end
        
        setParameters(transfo, trParams);
        
    end
end

end % classdef
