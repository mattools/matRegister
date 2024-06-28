classdef Transform < handle
%TRANSFORM  Abstract class for transform.
%   
%   Transform are designed to transform points into other points. This
%   class defines some abstract methods that have to be implemented by
%   derived classes.
%
%   Abstract classes:
%   transformPoint  - Computes coordinates of transformed point
%   transformVector - Computes coordinates of transformed vector
%   jacobianMatrix  - Computes jacobian matrix 
%
%   Example
%     % apply a tranlation transform to the point [10 10]
%     transfo = Translation([2 3]);
%     transformPoint(transfo, [10 10])
%     ans =
%         12    13
%
%   See also
%     AffineTransform, ParametricTransform

% ------
% Author: David Legland
% e-mail: david.legland@inra.fr
% Created: 2010-04-09,    using Matlab 7.9.0.529 (R2009b)
% Copyright 2010 INRA - Cepia Software Platform.

properties
    % the set of inner parameters of the transform
end


%% Abstract methods
methods (Abstract)
    
    % Return the dimension of this transform
    % In case of a projection transform, returns the dimension of input
    % points.
    getDimension(obj)
    
    % Computes coordinates of transformed point
    % pt2 = transformPoint(transfo, pt);
    transformPoint(obj, point)
        
    % Computes jacobian matrix, i.e. derivatives wrt to each coordinate
    % jacob(i,j) = d x_i / d x_j
    jacMat = jacobianMatrix(obj, position)
       
end % abstract methods

methods
    function jacMat = getJacobian(obj, position)
        % deprecated: use jacobianMatrix instead
        warning('deprecated: use method jacobianMatrix instead');
        jacMat = jacobianMatrix(obj, position);
    end
end

%% General methods
methods
    function res = transformVector(obj, vector, position)
        % Computes coordinates of transformed vector at a given position.
        %
        % vec2 = transformVector(transfo, vec, pos);
        %
        
        jac = jacobianMatrix(obj, position); % 2-by-2-by-N
        nv = size(vector, 1);
        res = zeros(nv, 2);
        for i = 1:nv
            res(i, :) = vector(i,:) * jac(:,:,i)';
        end
    end
    
    function res = compose(obj, that)
        % Computes the composition of the two transforms
        %
        % The following:
        % T = T1.compose(T2)
        % P2 = T.tansformPoint(P);
        % % is the same as 
        % P2 = T1.transformPoint(T2.transformPoint(P));
        res = ComposedTransform(that, obj);
    end
    
end % methods

%% Serialization methods
methods
    function write(obj, fileName, varargin)
        % Writes transform into a JSON file.
        % 
        % Requires implementation of the "toStruct" method.
        
        if exist('savejson', 'file') == 0
            error('Requires the ''jsonlab'' library');
        end
        if ~ismethod(obj, 'toStruct')
            error('Requires implementation of the ''toStruct'' method');
        end
        
        savejson('', toStruct(obj), 'FileName', fileName, varargin{:});
    end
end

methods (Static)
    function transfo = fromStruct(str)
        % Creates a new transform instance from a structure.
        
        % check existence of 'type' field
        if isfield(str, 'Type')
            type = str.Type;
        elseif isfield(str, 'type')
            type = str.type;
        else
            error('Requires a field with name "type"');
        end

        % parse transform
        try
            transfo = eval([type '.fromStruct(str)']);
        catch ME
            error(['Unable to parse transform with type: ' type]);
        end
    end
    
    function transfo = read(fileName)
        % Reads a transform from a file in JSON format.
        if exist('loadjson', 'file') == 0
            error('Requires the ''jsonlab'' library');
        end
        transfo = Transform.fromStruct(loadjson(fileName));
    end
end

end% classdef
