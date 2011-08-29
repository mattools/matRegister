classdef OptimizationIterationLogger < OptimizationListener
%OPTIMIZATIONITERATIONLOGGER  One-line description here, please.
%
%   Class OptimizationIterationLogger
%
%   Example
%   OptimizationIterationLogger
%
%   See also
%
%
% ------
% Author: David Legland
% e-mail: david.legland@grignon.inra.fr
% Created: 2011-08-29,    using Matlab 7.9.0.529 (R2009b)
% Copyright 2011 INRA - Cepia Software Platform.


%% Properties
properties
    % the file to write in (file handle)
    file;
    
    % time spent last iteration
    time0;
end % end properties


%% Constructor
methods
    function this = OptimizationIterationLogger(fileName)
        % Constructor for OptimizationIterationLogger class
        
        if nargin < 1
            error('Requires a file name as first argument');
        end
        
        this.file = fopen(fileName, 'wt');
    end

end % end constructors


%% Optimization Listener Methods
methods
    function optimizationStarted(this, src, event) %#ok<*INUSD>
        fprintf(this.file, 'iter     f_value     ||gradient||   time[s]\n');
        this.time0 = tic;
    end
    
    function optimizationIterated(this, optim, event)
        value = optim.value;
        time = toc(this.time0);
        grad = sqrt(sum(optim.gradient .^ 2));
        
        pattern = '%3d    %f    %f    %f\n';
        fprintf(this.file, pattern, optim.iter, value, grad, time);
        
        this.time0 = tic;
    end
    
    function optimizationTerminated(this, src, event)
        disp('End of optimization');
        fclose(this.file);
    end
end % end methods

end % end classdef

