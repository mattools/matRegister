classdef MatrixAffineTransform < AffineTransform
%MATRIXAFFINETRANSFORM  An affine transform defined by its matrix.
%
%   Example
%     % Creates a 2D transform corresponding to uniform scaling in x and y
%     mat = [2 0 0;0 2 0;0 0 1];
%     transfo = MatrixAffineTransform(mat);
%     transformPoint(transfo, [1 2])
%     ans =
%          2     4
%
%   See also
%     AffineTransform

% ------
% Author: David Legland
% e-mail: david.legland@inra.fr
% Created: 2010-06-17,    using Matlab 7.9.0.529 (R2009b)
% Copyright 2010 INRA - Cepia Software Platform.


%% Declaration of class properties
properties
    % center of the transform
    Matrix = eye(3);
end

%% Constructors
methods
    function obj = MatrixAffineTransform(varargin)
        % Create a new model for translation transform model
        
        if isempty(varargin)
            % parameters already set to default values
            return;
        end
        
        % extract first argument, and try to interpret
        var = varargin{1};
        if isa(var, 'MatrixAffineTransform')
            % copy constructor
            obj.Matrix = var.Matrix;
            
        elseif isnumeric(var)
            % initialisation constructor
            dim = size(var);
            if dim(1) == dim(2)
                obj.Matrix = var;
            elseif dim(1) == dim(2)-1
                % make the matrix square
                var(end+1, end) = 1;
                obj.Matrix = var;
            else
                error('Input matrix should be square');
            end
        else
            error('Unable to understand input arguments');
        end
    end % constructor declaration
end % methods

%% General methods
methods
    function mat = affineMatrix(obj)
        % Simply returns the inner matrix stored by obj transform.
        %
        mat = obj.Matrix;
    end
    
end % methods


%% Serialization methods
methods
    function str = toStruct(obj)
        % Converts to a structure to facilitate serialization
        str = struct('type', 'MatrixAffineTransform', ...
            'Matrix', obj.Matrix);
    end
end
methods (Static)
    function transfo = fromStruct(str)
        % Creates a new instance from a structure
        transfo = MatrixAffineTransform(str.Matrix);
    end
end

    
end % classdef