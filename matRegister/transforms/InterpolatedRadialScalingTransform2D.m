classdef InterpolatedRadialScalingTransform2D < Transform
%INTERPOLATEDRADIALSCALINGTRANSFORM2D Radial scaling transform in 2D.
%
%   Class InterpolatedRadialScalingTransform2D
%
%   Example
%   InterpolatedRadialScalingTransform2D
%
%   See also
%     RadialScalingTransform2D

% ------
% Author: David Legland
% e-mail: david.legland@inra.fr
% Created: 2018-12-04,    using Matlab 8.6.0.267246 (R2015b)
% Copyright 2018 INRA - BIA-BIBS.


%% Properties
properties
    % the center 
    Center = [0 0];
    
    % the angles for which sclaing is defined, in degrees, as 1-by-N array
    % default: 0:359
    Angles;
    
    % the scaling factor for each angle
    Scalings;
    
end % end properties


%% Constructor
methods
    function obj = InterpolatedRadialScalingTransform2D(varargin)
        % Constructor for InterpolatedRadialScalingTransform2D class
        %
        % T = InterpolatedRadialScalingTransform2D();
        % Empty constructor, using 360 scaling factors by default.
        %
        % T = InterpolatedRadialScalingTransform2D(N);
        % Creates a new transform by specifying the grid parameters.
        %
        % T = InterpolatedRadialScalingTransform2D(ANGLES);
        %
        % T = InterpolatedRadialScalingTransform2D(ANGLES, SCALINGS);
        
        if nargin == 0
            obj.Angles = 0:359;
                
        elseif nargin == 1
            var1 = varargin{1};
            if isa(var1, 'InterpolatedRadialScalingTransform2D')
                % copy constructor
                obj.Center = var1.Center;
                obj.Angles = var1.Angles;
                obj.Scalings = var1.Scalings;

            elseif isscalar(var1)
                n = var1;
                obj.Angles = linspace(0, 360, n+1);
                obj.Angles = obj.Angles(1:end-1);
            else
                obj.Angles = var1;
            end
            
        elseif nargin == 2
            obj.Angles = varargin{1};
            obj.Scalings = varargin{2};
        end

        if isempty(obj.Scalings)
            obj.Scalings = ones(size(obj.Angles));
        end
    end

end % end constructors


%% Methods implementing the Transform interface
methods
    function point2 = transformPoint(obj, point)
        % Compute coordinates of transformed point
        
        % polar coordinates wrt center point
        [theta, rho] = cart2pol(point(:,1)-obj.Center(1), point(:,2)-obj.Center(2));
        
        % convert to degrees
        ang = mod(rad2deg(theta) + 360, 360);
        
        % for each point, compute index of angle just lower
        inds0 = zeros(size(ang));
        for i = 1:length(ang)
            inds0(i) = find(ang(i) >= obj.Angles, 1, 'last');
        end
        
        % index of ange just after
        inds1 = mod(inds0, length(obj.Angles)) + 1;
        
        % angular position within angular interval
        t = (ang - obj.Angles(inds0)') / (obj.Angles(2) - obj.Angles(1));
        
        scals = (1-t) .* (obj.Scalings(inds0))' + t .* (obj.Scalings(inds1))';
        
        % apply scaling
        rho2 = rho .* scals;
        
        % convert back to cartesian coordinates
        [x2, y2] = pol2cart(theta, rho2);
        point2 = [x2+obj.Center(1) y2+obj.Center(2)];
    end
    
    function transformVector(obj, varargin)
        error('MatRegister:UnimplementedMethod', ...
            'Method "%s" is not implemented for class "%s"', ...
            'transformVector', mfilename);
    end
    
    function jacobianMatrix(obj, point) %#ok<INUSD>
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

    function dim = getDimension(obj) %#ok<MANU>
        dim = 2;
    end
end


%% Serialization methods
methods
    function str = toStruct(obj)
        % Converts to a structure to facilitate serialization
        str = struct('Type', 'InterpolatedRadialScalingTransform2D', ...
            'Angles', obj.Angles, ...
            'Scalings', obj.Scalings);
        if sum(obj.Center ~= 0) > 0
            str.Center = obj.Center;
        end
    end
end
methods (Static)
    function transfo = fromStruct(str)
        % Creates a new instance from a structure
        transfo = InterpolatedRadialScalingTransform2D(str.Angles, str.Scalings);
        if isfield(str, 'Center')
            transfo.Center = str.Center;
        end
    end
end

end % end classdef


