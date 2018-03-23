classdef BasicPolygonToPolygonMetric < handle
%BASICPOLYGONTOPOLYGONMETRIC  basic implementation of polygon to polygon comparator
%
%   output = BasicPolygonToPolygonMetric(input)
%
%   Example
%   BasicPolygonToPolygonMetric
%
%   See also
%
 
% ------
% Author: David Legland
% e-mail: david.legland@inra.fr
% Created: 2018-03-23,    using Matlab 9.3.0.713579 (R2017b)
% Copyright 2018 INRA - Cepia Software Platform.

properties
    % the reference polygon
    poly1;
    
    % the polygon to transform
    poly2;
    
    % the parametric transform
    transfo;
end

%% Constructor method
methods
    
    function this = BasicPolygonToPolygonMetric(poly1, poly2, transfo)
        
        if nargin < 3
            error('Require three input arguments');
        end
        
        if  size(poly1,2) ~= 2 || size(poly2,2) ~= 2
            error('Require 2D polygons as inputs');
        end
        
        % initialize fields
        this.poly1 = poly1;
        this.poly2 = poly2;
        this.transfo = transfo;
        
    end % constructor
end

%% Evaluation method
methods
    function that = evaluate(this, params)
        
        % update transformation model
        this.transfo.params = params;
        
        % transform polygon
        polyt = transformPoint(this.transfo, this.poly2);
        
        % compute metric
        dist = distancePointPolygon(this.poly1, polyt);
        that = sum(dist.^2);
    end
    
end
    
end
