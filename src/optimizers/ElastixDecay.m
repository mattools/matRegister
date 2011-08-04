classdef ElastixDecay < DecayFunction
%ELASTIXDECAY Decay function used in Elastix software
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
    a_num;
    a_denom = 50;
    alpha = 0.602;
end % end properties


%% Constructor
methods
    function this = ElastixDecay(varargin)
        % Constructor for ElastixDecay class
        
        if nargin > 0
            var = varargin{1};
            if isa(var, 'ElastixDecay')
                % Copy constructor
                this.a_num      = var.a_num;
                this.a_denom    = var.a_denom;
                this.alpha      = var.alpha;
                
            elseif isnumeric(var)
                % initialisation constructor
                this.a_num = var;
                
            else
                error('Wrong type of input');
                
            end
            varargin(1) = [];
            
            if ~isempty(varargin)
                this.a_denom = varargin{1};
            end
            if length(varargin) > 1
                this.alpha = varargin{2};
            end
        end
    end

end % end constructors


%% Methods
methods
    function alpha = evaluate(this, iter)
        alpha = this.a_num ./ (this.a_denom + iter) .^ this.alpha;
    end
end % end methods

end % end classdef

