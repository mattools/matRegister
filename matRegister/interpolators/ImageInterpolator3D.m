classdef ImageInterpolator3D < ImageInterpolator
% Abstract class that groups 3D image interpolators.
%   output = ImageInterpolator3D(input)
%

% ------
% Author: David Legland
% e-mail: david.legland@inra.fr
% Created: 2010-04-08,    using Matlab 7.9.0.529 (R2009b)
% Copyright 2010 INRA - Cepia Software Platform.

methods (Access = protected)
       function obj = ImageInterpolator3D(image)
        % Constructs a new ImageInterpolator3D object.
        
        % call superclass constructor
        obj = obj@ImageInterpolator(image);

    end % constructor declaration    
end
    
end % classdef

