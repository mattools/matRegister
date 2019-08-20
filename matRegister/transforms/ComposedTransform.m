classdef ComposedTransform < Transform
% Compose several transforms to create a new transform.
%   output = ComposedTransform(input)
%
%   Example
%   ParametricTransform2D
%
%   See also
%

% ------
% Author: David Legland
% e-mail: david.legland@inra.fr
% Created: 2010-04-09,    using Matlab 7.9.0.529 (R2009b)
% Copyright 2010 INRA - Cepia Software Platform.

properties
    % the set of encapsulated transforms
    Transforms;
end

%% Constructor
methods
    function obj = ComposedTransform(varargin)
        % class constructor
        % T = ComposedTransform(T1, T2);
        
        if isa(varargin{1}, 'ComposedTransform')
            % copy constructor
            var = varargin{1};
            obj.Transforms = var.Transforms;
        elseif isa(varargin{1}, 'Transform')
            % initialize tranform array
            nTransfos = length(varargin);
            obj.Transforms = cell(nTransfos, 1);
            for i = 1:nTransfos
                obj.Transforms{i} = varargin{i};
            end
        else
            error('Wrong parameter when constructing a Composed transform');
        end
    end % constructor
end


%% Methods implementing Transform interface
methods
    function dim = getDimension(obj)
        dim = getDimension(obj.Transforms{1});
    end

    function point = transformPoint(obj, point)
        for i = 1:length(obj.Transforms)
            point = transformPoint(obj.Transforms{i}, point);
        end
    end
    
    function vector = transformVector(obj, vector, position)
        for i = 1:length(obj.Transforms)
            vector = transformVector(obj.Transforms{i}, vector, position);
        end
    end
    
    function jacMat = jacobianMatrix(obj, point)
        % Compute jacobian matrix, i.e. derivatives for coordinate
        % jacob(i,j) = d x_i / d x_j
        jacMat = jacobianMatrix(obj.Transforms{1}, point);
        for i = 2:length(obj.Transforms)
            jacMat = jacobianMatrix(obj.Transforms{i}, point) * jacMat;
        end
    end
end % methods

end % classdef
