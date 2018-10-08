classdef Transform < handle
%TRANSFORM  Abstract class for transform
%   
%   Transform are designed to transform points into other points. This
%   class defines some abstract methods that have to be implemented by
%   derived classes.
%
%   Abstract classes:
%   transformPoint  - Computes coordinates of transformed point
%   transformVector - Computes coordinates of transformed vector
%   getJacobian     - Computes jacobian matrix 
%
%   Example
%   trans = (...); % define a transform by using a derived class 
%   pt = (...);    % create a point corresponding to transform input
%   pt2 = trans.transformPoint(pt);
%
%   See also
%
%
% ------
% Author: David Legland
% e-mail: david.legland@grignon.inra.fr
% Created: 2010-04-09,    using Matlab 7.9.0.529 (R2009b)
% Copyright 2010 INRA - Cepia Software Platform.

properties
    % the set of inner parameters of the transform
end

%% Static methods
methods (Static)
end

%% Abstract methods
methods (Abstract)
    
    getDimension(this)
    % GETDIMENSION Return the dimension of this transform
    % In case of a projection transform, returns the dimension of input
    % points.
    
    transformPoint(this, point)
    % TRANSFORMPOINT Computes coordinates of transformed point
    % PT2 = this.transformPoint(PT);
    
    transformVector(this, vector, position)
    % TRANSFORMVECTOR Computes coordinates of transformed vector
    % VEC2 = this.transformPoint(VEC, PT);
    
    jacMat = jacobianMatrix(this, position)
    % Computes jacobian matrix, i.e. derivatives wrt to each coordinate
    % jacob(i,j) = d x_i / d x_j
       
end % abstract methods

methods
    function jacMat = getJacobian(this, position)
        warning('deprecated: use method jacobianMatrix instead');
        jacMat = jacobianMatrix(this, position);
    end
end

%% General methods
methods
    
    function res = compose(this, that)
        % Computes the composition of the two transforms
        %
        % The following:
        % T = T1.compose(T2)
        % P2 = T.tansformPoint(P);
        % % is the same as 
        % P2 = T1.transformPoint(T2.transformPoint(P));
        res = ComposedTransform(that, this);
    end
    
end % methods

%% Serialization methods
methods
    function write(this, fileName, varargin)
        % Writes transform into a JSON file
        % 
        % Requires implementation of the "toStruct" method.
        
        if exist('savejson', 'file') == 0
            error('Requires the ''jsonlab'' library');
        end
        savejson('', toStruct(this), 'FileName', fileName, varargin{:});
    end
end

methods (Static)
    function transfo = fromStruct(str)
        % Creates a new transform instance from a structure
        
        % check existence of 'type' field
        if ~isfield(str, 'type')
            error('Requires a field with name "type"');
        end
        type = str.type;

        % parse transform
        try
            transfo = eval([type '.fromStruct(str)']);
        catch ME
            error(['Unable to parse transform with type: ' type]);
        end
    end
    
    function transfo = read(fileName)
        % Reads a transform from a file in JSON format
        if exist('loadjson', 'file') == 0
            error('Requires the ''jsonlab'' library');
        end
        transfo = Transform.fromStruct(loadjson(fileName));
    end
end

end% classdef
