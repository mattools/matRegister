classdef RadialScalingTransform2D < Transform
%RADIALSCALINGTRANSFORM2D Radial scaling transform in 2D
%
%   Class RadialScalingTransform2D
%
%   Example
%   RadialScalingTransform2D
%
%   See also
%

% ------
% Author: David Legland
% e-mail: david.legland@inra.fr
% Created: 2018-11-11,    using Matlab 8.6.0.267246 (R2015b)
% Copyright 2018 INRA - BIA-BIBS.


%% Properties
properties
    % the center 
    center = [0 0];
    
    % the angles for which sclaing is defined, in degrees, as 1-by-N array
    % default: 0:359
    angles;
    
    % the scaling factor for each angle
    scalings;
    
end % end properties


%% Constructor
methods
    function this = RadialScalingTransform2D(varargin)
        % Constructor for RadialScalingTransform2D class
        %
        % T = RadialScalingTransform2D();
        % Empty constructor, using 360 scaling factors by default.
        %
        % T = RadialScalingTransform2D(N);
        % Creates a new transform by specifying the grid parameters.
        %
        % T = RadialScalingTransform2D(ANGLES);
        %
        % T = RadialScalingTransform2D(ANGLES, SCALINGS);
        
        if nargin == 0
            this.angles = 0:359;
                
        elseif nargin == 1
            var1 = varargin{1};
            if isa(var1, 'RadialScalingTransform2D')
                % copy constructor
                this.center = var1.center;
                this.angles = var1.angles;
                this.scalings = var1.scalings;

            elseif isscalar(var1)
                n = var1;
                this.angles = linspace(0, 360, n+1);
                this.angles = this.angles(1:end-1);
            else
                this.angles = var1;
            end
            
        elseif nargin == 2
            this.angles = varargin{1};
            this.scalings = varargin{2};
        end

        if isempty(this.scalings)
            this.scalings = ones(size(this.angles));
        end
    end

end % end constructors


%% Methods implementing the Transform interface
methods
    function point2 = transformPoint(this, point)
        % Compute coordinates of transformed point
        
        [theta, rho] = cart2pol(point(:,1)-this.center(1), point(:,2)-this.center(2));
        ang = mod(rad2deg(theta) + 360, 360);
        
        inds = zeros(size(ang));
        for i = 1:length(ang)
            [tmp, inds(i)] = min((ang(i) - this.angles).^2); %#ok<ASGLU>
        end
        
        rho2 = rho .* (this.scalings(inds))';
        
        [x2, y2] = pol2cart(theta, rho2);
        point2 = [x2+this.center(1) y2+this.center(2)];
    end
    
    function transformVector(this, varargin)
        error('MatRegister:UnimplementedMethod', ...
            'Method "%s" is not implemented for class "%s"', ...
            'transformVector', mfilename);
    end
    
    function jacobianMatrix(this, point) %#ok<INUSD>
        % Jacobian matrix of the given point
        %
        %   JAC = jacobianMatrix(TRANS, PT)
        %   where PT is a N-by-2 array of points, returns the spatial
        %   jacobian matrix of each point in the form of a 2-by-2-by-N
        %   array.
        %
        
        error('MatRegister:UnsupportedOperation', ...
            'Method "%s" is not supported for class "%s"', ...
            'getJacobian', mfilename);
    end

    function dim = getDimension(this) %#ok<MANU>
        dim = 2;
    end
end


%% Serialization methods
methods
    function str = toStruct(this)
        % Converts to a structure to facilitate serialization
        str = struct('type', 'RadialScalingTransform2D', ...
            'angles', this.angles, ...
            'scalings', this.scalings);
        if sum(this.center ~= 0) > 0
            str.center = this.center;
        end
    end
end
methods (Static)
    function transfo = fromStruct(str)
        % Creates a new instance from a structure
        transfo = RadialScalingTransform2D(str.angles, str.scalings);
        if isfield(str, 'center')
            transfo.center = str.center;
        end
    end
end

end % end classdef

