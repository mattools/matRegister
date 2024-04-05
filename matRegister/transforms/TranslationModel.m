classdef TranslationModel < ParametricTransform & AffineTransform
%Transformation model for a translation defined by ND parameters.
%   output = TranslationModel(input)
%
%   Parameters of the transform:
%   inner optimisable parameters of the transform have following form:
%   params(1): tx       (in user spatial unit)
%   params(2): ty       (in user spatial unit)
%
%

% ------
% Author: David Legland
% e-mail: david.legland@inra.fr
% Created: 2010-02-17,    using Matlab 7.9.0.529 (R2009b)
% Copyright 2010 INRA - Cepia Software Platform.


%% Constructors
methods
    function obj = TranslationModel(varargin)
        % Create a new model for translation transform model.
        
        % set parameters to default translation in 2D
        obj.Params = [0 0];
        
        if ~isempty(varargin)
            % extract first argument, and try to interpret
            var = varargin{1};
            if isa(var, 'TranslationModel')
                % copy constructor
                obj.Params = var.Params;
                
            elseif isnumeric(var)
                if isscalar(var)
                    % interpret the scalar as the working dimension
                    obj.Params = zeros(1, var);
                elseif size(var, 1)==1
                    % interpret the vector as the translation parameter
                    obj.Params = var;
                else
                    error('Input argument must be a scalar or a row vector');
                end
                
            else
                error('Unable to understand input arguments');
            end
        end
        
        % update parameter names
        np = length(obj.Params);
        switch np
            case 2
                obj.ParamNames = {'X shift', 'Y shift'};
            case 3
                obj.ParamNames = {'X shift', 'Y shift', 'Z shift'};
            otherwise
                obj.ParamNames = cellstr(num2str((1:4)', 'Shift %d'));
        end
        
    end % constructor declaration
end


%% Implementation of methods from ParametricTransform 
methods
    function jac = parametricJacobian(obj, x, varargin) %#ok<INUSD>
        % Compute jacobian matrix, i.e. derivatives for each parameter.
        nd = length(obj.Params);
        jac = eye(nd);
    end
    
    function transfo = clone(obj)
        transfo = TranslationModel(obj.Params);
    end

end % methods


%% Implementation of AffineTransform methods
methods
    function mat = affineMatrix(obj)
        % Returns the affine matrix that represents obj transform.
        nd = length(obj.Params);
        mat = eye(nd+1);
        mat(1:end-1, end) = obj.Params(:);
    end
end


%% Serialization methods
methods
    function str = toStruct(obj)
        % Converts to a structure to facilitate serialization.
        str = struct('Type', 'TranslationModel', ...
            'Parameters', obj.Params);
    end
end
methods (Static)
    function motion = fromStruct(str)
        % Creates a new instance from a structure.
        params = str.Parameters;
        motion = TranslationModel(params);
    end
end


end % classdef