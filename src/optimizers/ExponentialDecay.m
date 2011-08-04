classdef ExponentialDecay < DecayFunction
%EXPONENTIALDECAY Exponential decay function.
%
%   DEC = ExponentialDecay(TAU)
%   DEC = ExponentialDecay(TAU, ALPHA0)
%   Creates a new ExponentialDecay function, whose value is evaluated from:
%     dec(t) = alpha0 * exp(- t /tau)
%   
%   Example
%     dec = ExponentialDecay(200);
%     t = 1:500;
%     plot(t, evaluate(dec, t));
%     title('Exponential decay, tau=200');
%
%   See also
%     ElastixDecay, DecayFunction
%
%
% ------
% Author: David Legland
% e-mail: david.legland@grignon.inra.fr
% Created: 2011-08-04,    using Matlab 7.9.0.529 (R2009b)
% Copyright 2011 INRA - Cepia Software Platform.


%% Properties
properties
    tau = 1;
    alpha0 = 1;
end % end properties


%% Constructor
methods
    function this = ExponentialDecay(varargin)
        % Constructor for ExponentialDecay class
        
        if nargin > 0
            var = varargin{1};
            if isa(var, 'ExponentialDecay')
                this.tau = var.tau;
                
            elseif isnumeric(var)
                this.tau = var;
                
            else
                error('Wrong type of input');
                
            end
            varargin(1) = [];
            
            if ~isempty(varargin)
                this.alpha0 = varargin{1};
            end
        end
    end

end % end constructors


%% Methods
methods
    function alpha = evaluate(this, iter)
        alpha = this.alpha0 * exp(-iter ./ this.tau);
    end
end % end methods

end % end classdef

