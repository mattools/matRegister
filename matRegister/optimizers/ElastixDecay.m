classdef ElastixDecay < DecayFunction
% Decay function used in Elastix software.
%
%   Class ElastixDecay
%
%   a_k = a0 / (A + k)^alpha
%   (we drop the +1 from original Elastix convention, as indexing is
%   usually different with Matlab)
%
%   Example
%     dec = ElastixDecay(1000, 50, 0.602);
%     t = 1:500;
%     plot(t, evaluate(dec, t));
%     title('Elastix Decay, a=1000, A=50, alpha=0.602');
%
%
%   See also
%     ExponentialDecay, DecayFunction
%
% ------
% Author: David Legland
% e-mail: david.legland@grignon.inra.fr
% Created: 2011-08-04,    using Matlab 7.9.0.529 (R2009b)
% Copyright 2011 INRA - Cepia Software Platform.


%% Properties
properties
    A_num;
    A_denom = 50;
    Alpha = 0.602;
end % end properties


%% Constructor
methods
    function obj = ElastixDecay(varargin)
        % Constructor for ElastixDecay class
        
        if nargin > 0
            var = varargin{1};
            if isa(var, 'ElastixDecay')
                % Copy constructor
                obj.A_num      = var.A_num;
                obj.A_denom    = var.A_denom;
                obj.Alpha      = var.Alpha;
                
            elseif isnumeric(var)
                % initialisation constructor
                obj.A_num = var;
                
            else
                error('Wrong type of input');
                
            end
            varargin(1) = [];
            
            if ~isempty(varargin)
                obj.A_denom = varargin{1};
            end
            if length(varargin) > 1
                obj.Alpha = varargin{2};
            end
        end
    end

end % end constructors


%% Methods
methods
    function alpha = evaluate(obj, iter)
        alpha = obj.A_num ./ (obj.A_denom + iter) .^ obj.Alpha;
    end
end % end methods

end % end classdef

