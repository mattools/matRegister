classdef OptimizationIterationLogger < OptimizationListener
% Write iteration info into a log file.
%
%   Class OptimizationIterationLogger
%
%   Example
%   OptimizationIterationLogger
%
%   See also
%

% ------
% Author: David Legland
% e-mail: david.legland@inra.fr
% Created: 2011-08-29,    using Matlab 7.9.0.529 (R2009b)
% Copyright 2011 INRA - Cepia Software Platform.


%% Properties
properties
    % the file to write in (file handle)
    File;
    
    % time spent in last iteration
    Time0;
    
end % end properties


%% Constructor
methods
    function obj = OptimizationIterationLogger(fileName)
        % Constructor for OptimizationIterationLogger class
        
        if nargin < 1
            error('Requires a file name as first argument');
        end
        
        obj.File = fopen(fileName, 'wt');
    end

end % end constructors


%% Optimization Listener Methods
methods
    function optimizationStarted(obj, src, event) %#ok<*INUSD>
        fprintf(obj.File, 'iter     f_value     ||gradient||   time[s]\n');
        obj.Time0 = tic;
    end
    
    function optimizationIterated(obj, optim, event)
        value = optim.Value;
        time = toc(obj.Time0);
        grad = sqrt(sum(optim.Gradient .^ 2));
        
        pattern = '%3d    %f    %f    %f\n';
        fprintf(obj.File, pattern, optim.Iter, value, grad, time);
        
        obj.Time0 = tic;
    end
    
    function optimizationTerminated(obj, src, event)
        disp('End of optimization');
        fclose(obj.File);
    end
end % end methods

end % end classdef

