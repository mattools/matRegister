classdef BackwardTransformedImage < ImageFunction
% Describes an image associated with a backward transform.
%
%   TIM = BackwardTransformedImage(IMG, TRANS, 'linear');
%   where IMG is an image, and TRANS a transform.
%
%   TIM = BackwardTransformedImage(interpolator, TRANS);
%   where interpolator is an instance of ImageInterpolator
%
%   Example
%   img = Image2D('cameraman.tif');
%   interpolator = LinearInterpolator2D(img);
%   trans = Translation2D(30, 60);
%   tim = BackwardTransformedImage(interpolator, trans);
%   rsp = ImageResampler(img);
%   img2 = rsp.resample(tim);
%   img2.show();
%
%   See also
%
%

% ------
% Author: David Legland
% e-mail: david.legland@inra.fr
% Created: 2010-04-08,    using Matlab 7.9.0.529 (R2009b)
% Copyright 2010 INRA - Cepia Software Platform.
    
%% Declaration of class properties
properties        
    % the transform from physical space to image space
    Transform       = [];

    % an interpolator, that must inherit ImageFunction
    Interpolator    = [];
end

%% Static methods
methods(Static)
    function tim = create(varargin)
        tim = BackwardTransformedImage(varargin{:});
    end
end

%% Constructors
methods
    function obj = BackwardTransformedImage(varargin)
        % Constructs a new BackwardTransformedImage object.
        if nargin==0
            % empty constructor
            % (nothing to do !)

        elseif nargin==1 && isa(varargin{1}, 'BackwardTransformedImage')
            % copy constructor: copy each field
            tim = varargin{1};
            obj.Transform      = tim.Transform;
            obj.Interpolator   = tim.Interpolator;

        elseif isa(varargin{1}, 'ImageFunction')
            % standard constructor: interpolator and transform
            obj.Interpolator   = varargin{1};
            obj.Transform      = varargin{2};

        elseif isa(varargin{1}, 'Image')
            % can also specify image+transform+interpolator type
            img = varargin{1};
            obj.Transform = varargin{2};
            % parse the interpolator
            if length(varargin)>2
                var = varargin{3};
                if ischar(var)
                    obj.Interpolator = ImageInterpolator.create(img, var);
                elseif isa(var, 'ImageFunction')
                    obj.Interpolator = var;
                end
            else
                % if no interpolation method is specified, use linear
                % interpolator as default
                obj.Interpolator = ImageInterpolator.create(img, 'linear');
            end

        else
            error('Error in input arguments');
        end
    end % constructor declaration
end

end % classdef

