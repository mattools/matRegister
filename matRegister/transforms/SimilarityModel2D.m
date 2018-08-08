classdef SimilarityModel2D < AffineTransform & ParametricTransform
%Transformation model for a similarity: rotation+scaling followed by a translation
%   
%   Inner optimisable parameters of the transform have the following form:
%   params[1] = tx
%   params[2] = ty
%   params[3] = theta, in degrees
%   params[4] = logk, the binary logarithm of the scaling factor
%       ( 0 -> no scaling)
%       (+1 -> uniform scaling by a factor of 2)
%       (-1 -> uniform scaling by a factor of 1/2)
% 
%   TRANS = SimilarityModel2D();
%   Create using default parameters (all zero).
%
%   TRANS = SimilarityModel2D(PARAMS);
%   Create a new transform model by initializing parameters. PARAMS is a
%   1-by-4 row vector containing TX, TY, THETA and LOGK parameters.
%
%

% ------
% Author: David Legland
% e-mail: david.legland@inra.fr
% Created: 2017-08-07,    using Matlab 7.9.0.529 (R2009b)
% Copyright 2010 INRA - Cepia Software Platform.


%% Constructors
methods
    function this = SimilarityModel2D(varargin)
        % Create a new model for 2D motion transform model
        % (defined by 2 translation parameters and 1 rotation angle)

        % call parent constructor for initializing parameters
        this = this@ParametricTransform([0 0 0 0]);
        
        % setup default parameters
        this.params = [0 0 0 0];

        if ~isempty(varargin)
            % extract first argument, and try to interpret
            var = varargin{1};
            if isa(var, 'SimilarityModel2D')
                this.params = var.params;
                
            elseif isnumeric(var)
                if length(var) == 1 && var(1) ~= 2
                    % if argument is a scalar, it corresponds to dimension
                    error('SimilarityModel2D is defined only for 2 dimensions');
                    
                elseif length(var) == 4
                    % if argument is a row vector, it is the parameters
                    this.params = var(1:4);
                else
                    error('Please specify translation, rotation angle and scaling');
                end
            else
                error('Unable to understand input arguments');
            end
        end
       
        % setup parameter names
        this.paramNames = {'X shift', 'Y shift', 'Theta (°)', 'Log2 scale'};
                
    end % constructor declaration
end

%% Standard methods
methods        
    function dim = getDimension(this) %#ok<MANU>
        dim = 2;
    end

    function mat = getAffineMatrix(this)
        % Compute affine matrix associated with this transform
        
        % translation vector
        tx = this.params(1);
        ty = this.params(2);
        
        % rotation angle, converted to radians
        theta = this.params(3) * pi / 180;

        % scaling factor, converted to ususal scale
        k = power(2, this.params(4));
        
        % pre-computations
        cot = cos(theta);
        sit = sin(theta);

        % build matrix
        mat = [...
        	k*cot -k*sit tx ; ...
        	k*sit +k*cot ty ; ...
              0      0    1 ];

    end
  
    function jacobian = getParametricJacobian(this, x, varargin)
        % Compute jacobian matrix, i.e. derivatives for each parameter
       
        % extract coordinate of input point(s)
        if isempty(varargin)
            y = x(:,2);
            x = x(:,1);
        else
            y = varargin{1};
        end
        
        % convert angles to radians
        theta = this.params(3) * pi / 180;

        % scaling factor, converted to ususal scale
        k = power(2, this.params(4));
        % the partial derivative of 2^k wrt k
        dk = log(2)*this.params(k);
       
        % precompute angle functions
        cot = cos(theta);
        sit = sin(theta);
        
        % compute parametric jacobian
        jacobian = [...
            1 0 k*(-sit*x - cot*y)*pi/180 (cot*x-sit*y)*dk ;
            0 1 k*( cot*x - sit*y)*pi/180 (sit*x+cot*y)*dk];
    end

end % methods

% %% I/O Methods
% methods (Static)
%     function transfo = readFromFile(fileName)
%         % Read transform from the given file name.
%         % Returns a new instance of MotionModel2D.
%         %
%         % Example
%         %   TRANSFO = MotionModel2D.readFromFile('transfo.txt');
%         
%         map = readPropertyFile(fileName);
%         transfo = MotionModel2D.createFromPropertyMap(map);
%     end
%     
%     function transfo = createFromPropertyMap(map)
%         % Create a new transform from a set of properties
%         
%         transfo = MotionModel2D();
%         
%         nbParams = str2double(map('TransformParameterNumber'));
%         
%         trParams = map('TransformParameters');
%         trParams= cellfun(@str2double, regexp(trParams, '\s*', 'split'));
%         
%         if nbParams ~= length(trParams)
%             error('Wrong number of parameters');
%         end
%         
%         setParameters(transfo, trParams);
%         
%         center = map('TransformCenter');
%         center = cellfun(@str2double, regexp(center, '\s*', 'split'));
%         transfo.center = center;
%     end
% end

end % classdef