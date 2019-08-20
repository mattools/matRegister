classdef AffineTransformModel2D < AffineTransform & ParametricTransform 
%Transformation model for a 2D affine transform.
%
%   TRANSFO = AffineTransformModel2D()
%   Creates a new instanfce of affine parametric model in 2D, with 6
%   parameters.
%
%   Example
%   AffineTransformModel2D
%
%   See also
%     AffineTransform

% ------
% Author: David Legland
% e-mail: david.legland@inra.fr
% Created: 2017-09-29,    using Matlab 7.9.0.529 (R2009b)
% Copyright 2017 INRA - Cepia Software Platform.


%% Constructors
methods
    function obj = AffineTransformModel2D(varargin)
        % Create a new affine transform model
        
        % initialize to identity
        obj.Params = [1 0 0   0 1 0];
        
        if ~isempty(varargin)
            % extract first argument, and try to interpret
            var = varargin{1};
            if isa(var, 'AffineTransformModel2D')
                % copy constructor
                obj.Params = var.Params;
                
            elseif isnumeric(var)
                len = length(var);
                if len == 6
                    % all parameters are specified as a row vector
                    obj.Params = var;
                    
                elseif sum(size(var)~=[3 3])==0 || sum(size(var)~=[2 3])==0
                    % the parameter matrix is specified
                    var = var(1:2, :)';
                    obj.Params = var(:)';
                    
                elseif len == 1
                    % give the possibility to call constructor with the
                    % dimension of image
                    % Initialize with identity transform
                    if var ~= 2
                        error('Defined only for 2 dimensions');
                    end
                    
                else
                    error('Please specify affine parameters');
                end
                
            else
                error('Unable to understand input arguments');
            end
        end
        
        % setup parameter names
        obj.ParamNames = {...
            'm00', 'm01', 'm02', ...
            'm10', 'm11', 'm12', ...
            };
    end     
    
end

%% Methods specific to obj class
methods
    function initFromAffineTransform(obj, transform)
        % Initialize parameters from an affine transform (class or matrix)
        %
        % Example
        %   T = AffineTransformModel2D();
        %   mat = createRotation([10 20], .3);
        %   T.initFromAffineTransform(mat);
        %   T.getParameters()
        %
        
        % if the first argument is a Transform, extract its affine matrix
        if isa(transform, 'AffineTransform')
            transform = affineMatrix(transform);
        end
        
        % format matrix to have a row vector of 12 elements
        params0 = transform(1:3, :)';
        obj.Params = params0(:)';
    end

end


%% Implementation of methods inherited from AffineTransform
methods
    function dim = getDimension(obj) %#ok<MANU>
        % Returns the dimension of the tansform, in obj case 2.
        dim = 2;
    end

    function mat = affineMatrix(obj)
        % Compute affine matrix associated with obj transform
        
        % convert parameters to a 3-by-3 affine matrix
        mat = [obj.Params(1:3); obj.Params(4:6); 0 0 1];
    end
end


%% Implementation of methods inherited from ParametricTransform
methods
    function jacobian = parametricJacobian(obj, x, varargin) %#ok<INUSL>
        % Compute jacobian matrix, i.e. derivatives for each parameter
        % Result is a 2-by-6 matrix
        
        % extract coordinate of input point(s)
        if isempty(varargin)
            y = x(:,2);
            x = x(:,1);
        else
            y = varargin{1};
        end
        
        % compute the Jacobian matrix using pre-computed elements
        jacobian = [...
            x y 1       zeros(1,3); ...
            zeros(1,3)  x y 1; ...
            ];
        
    end
    
end % parametric transform methods 


%% Serialization methods
methods
    function str = toStruct(obj)
        % Converts to a structure to facilitate serialization
        str = struct('Type', 'AffineTransformModel2D', ...
            'Parameters', obj.Params);
    end
end
methods (Static)
    function motion = fromStruct(str)
        % Creates a new instance from a structure
        params = str.Parameters;
        motion = AffineTransformModel2D(params);
    end
end


end
