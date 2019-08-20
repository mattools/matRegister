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
    Poly1;
    
    % the polygon to transform
    Poly2;
    
    % the parametric transform
    Transfo;
end

%% Constructor method
methods
    
    function obj = BasicPolygonToPolygonMetric(poly1, poly2, transfo)
        
        if nargin < 3
            error('Require three input arguments');
        end
        
        if  size(poly1,2) ~= 2 || size(poly2,2) ~= 2
            error('Require 2D polygons as inputs');
        end
        
        % initialize fields
        obj.Poly1 = poly1;
        obj.Poly2 = poly2;
        obj.Transfo = transfo;
        
    end % constructor
end

%% Evaluation method
methods
    function that = evaluate(obj, params)
        
        % update transformation model
        obj.Transfo.Params = params;
        
        % transform polygon
        polyt = transformPoint(obj.Transfo, obj.Poly2);
        
        % compute metric
        dist = distancePointPolygon(obj.Poly1, polyt);
        that = sum(dist.^2);
    end
    
end
    
end
