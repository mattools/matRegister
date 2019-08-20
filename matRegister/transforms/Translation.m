classdef Translation < AffineTransform
%TRANSLATION Defines a translation in ND space.
%   
%   Representation of a translation transform in ND space.
%   The translation vector is stored by the class, and the transformation
%   matrix is computed when needed.
%
%   For a parameterized translation, see TranslationModel class.
%
%   Example
%     % Apply 2D translation on cameraman image   
%     img = imread('cameraman.tif');
%     transfo = Translation([50 30]);
%     [x, y] = meshgrid(1:256, 1:256);
%     pts2 = transfo.transformPoint([x(:) y(:)]);
%     res = zeros(size(img), 'uint8');
%     res(:) = imEvaluate(img, pts2);
%     figure; imshow(res);

%   % Creates a 3D translation 
%     T = Translation([3 4 5]);
%
%   See also
%     AffineTransform, TranslationModel
%

% ------
% Author: David Legland
% e-mail: david.legland@inra.fr
% Created: 2010-04-09,    using Matlab 7.9.0.529 (R2009b)
% Copyright 2010 INRA - Cepia Software Platform.

properties
    % the translation vector, stored as a row vector
    U;
end

%% Abstract methods
methods
    function obj = Translation(varargin)
        %TRANSLATION Class constructor
        %
        % T = Translation(VEC)
        % where VEC is a row vector, initialize translation with vector VEC
        %
        % T = Translation(V1, V2...)
        % Initialize each component separately. Each component must be a
        % scalar.
        %
        % T = Translation;
        % Initialize with 2D translation with (0,0) vector.
        %
        % T = Translation(T0)
        % copy constructor from another Translation object/
        
        if nargin == 0
            % empty constructor, initialized to (0,0) 2D translation
            obj.U = [0 0];
            return;
        end
        
        if isa(varargin{1}, 'Translation')
            % copy constructor
            var = varargin{1};
            obj.U = var.U;
        elseif isnumeric(varargin{1})
            var = varargin{1};
            if size(var, 2) == 1
                % initialize with separate scalar parameters
                obj.U = zeros(1, nargin);
                for i=1:nargin
                    obj.U(i) = varargin{i};
                end
            else
                % initialize with bundled parameters
                obj.U = var;
            end
        else
            error('Wrong parameter when constructing a Translation');
        end
    end
end

methods
    function d = getDimension(obj)
        d = length(obj.U);
    end
    
    function mat = affineMatrix(obj)
        % Returns the (ND+1)*(ND+1) affine matrix representing translation
        nd = length(obj.U);
        mat = eye(nd+1);
        mat(1:end-1, end) = obj.U(:);
    end
    
end % end methods

%% Serialization methods
methods
    function str = toStruct(obj)
        % Converts to a structure to facilitate serialization
        str = struct('Type', 'Translation', ...
            'Vector', obj.U);
    end
end
methods (Static)
    function transfo = fromStruct(str)
        % Creates a new instance from a structure
        transfo = Translation(str.Vector);
    end
end

end % end classdef
