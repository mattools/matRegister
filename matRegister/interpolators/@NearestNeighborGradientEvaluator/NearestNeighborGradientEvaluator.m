classdef NearestNeighborGradientEvaluator < ImageFunction
%NearestNeighborGradientEvaluator  One-line description here, please.
%
%   output = NearestNeighborGradientEvaluator(input)
%
%   Example
%   NearestNeighborGradientEvaluator
%
%   See also
%
%

% ------
% Author: David Legland
% e-mail: david.legland@inra.fr
% Created: 2010-12-09,    using Matlab 7.9.0.529 (R2009b)
% Copyright 2010 INRA - Cepia Software Platform.

properties
    % the image whom gradient will be evaluated
    Image;
    
    % output type of interpolation, double by default
    OutputType = 'double';

    % default value for points outside image
    DefaultValue = NaN;
    
    % the filters used for gradient evaluation in each direction.
    % should be a cell array with as many columns as image dimension.
    Filters;
end

%% Constructors
methods
    function obj = NearestNeighborGradientEvaluator(varargin)
       
        if nargin == 0
            return;
        end
        
        % extract image to interpolate
        var = varargin{1};
        if isa(var, 'NearestNeighborGradientEvaluator')
            obj.Image = var.Image;
        elseif isa(var, 'Image')
            obj.Image = var;
        end
        varargin(1) = [];
        
        if channelNumber(obj.Image) > 1
            error('Gradient of vector images is not supported');
        end
        
        % setup default values depending on image dimension
        nd = ndims(obj.Image);
        switch nd
            case 1
                obj.Filters = {[1 0 -1]'};
            case 2
                sx = fspecial('sobel')'/8;
                obj.Filters = {sx, sx'};
            case 3
                [sx, sy, sz] = Image.create3dGradientKernels();
                obj.Filters = {sx, sy, sz};
        end
        
        % parse user specified options
        while length(varargin) > 1
            paramName = varargin{1};
            if strcmpi(paramName, 'defaultValue')
                obj.DefaultValue = varargin{2};
                
            elseif strcmpi(paramName, 'filters')
                obj.Filters = varargin{2};
                
            elseif strcmpi(paramName, 'outputType')
                obj.OutputType = varargin{2};
                
            else
                error(['Unknown parameter name: ' paramName]);
            end
            
            varargin(1:2) = [];
        end
    end
end

end