classdef MatrixAffineTransform < AffineTransform
%MATRIXAFFINETRANSFORM  An affine transform defined by its matrix
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
    matrix = eye(3);
end

%% Constructors
methods
    function this = MatrixAffineTransform(varargin)
        % Create a new model for translation transform model
        
        if isempty(varargin)
            % parameters already set to default values
            return;
        end
        
        % extract first argument, and try to interpret
        var = varargin{1};
        if isa(var, 'MatrixAffineTransform')
            % copy constructor
            this.matrix = var.matrix;
            
        elseif isnumeric(var)
            % initialisation constructor
            dim = size(var);
            if dim(1) == dim(2)
                this.matrix = var;
            elseif dim(1) == dim(2)-1
                % make the matrix square
                var(end+1, end) = 1;
                this.matrix = var;
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
    function mat = affineMatrix(this)
        % Simply returns the inner matrix stored by this transform.
        %
        mat = this.matrix;
    end
    
end % methods


%% Serialization methods
methods
    function str = toStruct(this)
        % Converts to a structure to facilitate serialization
        str = struct('type', 'MatrixAffineTransform', ...
            'matrix', this.matrix);
    end
end
methods (Static)
    function transfo = fromStruct(str)
        % Creates a new instance from a structure
        transfo = MatrixAffineTransform(str.matrix);
    end
end

    
end % classdef