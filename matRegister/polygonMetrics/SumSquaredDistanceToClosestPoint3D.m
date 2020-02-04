classdef SumSquaredDistanceToClosestPoint3D < handle
% SSD metric for two 3D point clouds.
%
%   output = SumSqDistancesVertexToClosestVertex3d(PTS1, PTS2, TRANSFO)
%
%   This version computes distances by iterating over the reference points,
%   and computing distance to closest point among transformed point set 2.
%
%   Example
%   SumSquaredDistanceToClosestPoint3D
%
%   See also
%
 
% ------
% Author: David Legland
% e-mail: david.legland@inrae.fr
% INRAE - BIA Research Unit - BIBS Platform (Nantes)
% Created: 2020-01-22,    using Matlab 9.7.0.1247435 (R2019b) Update 2
% Copyright 2020 INRAE.


properties
    % the reference points
    RefPoints;
    
    % the points to transform
    Points;
    
    % the parametric transform
    Transfo;
end

%% Constructor method
methods
    function obj = SumSquaredDistanceToClosestPoint3D(refPoints, points, transfo)
        % Creates a new metric based on vertex to vertex distances
        %
        % Input arguments are the reference vertex set, the moving set of
        % points, and the transform model.
        %
        % metric = SumDistancesVertexToClosestVertex3d2(refVerts, verts, transfo);
        
        if nargin < 3
            error('Require three input arguments');
        end
        
        if  size(refPoints, 2) ~= 3 || size(points, 2) ~= 3
            error('Require 3D points as inputs');
        end
        
        % initialize fields
        obj.RefPoints = refPoints;
        obj.Points = points;
        obj.Transfo = transfo;
        
    end % constructor
end

%% Evaluation method
methods
    function res = evaluate(obj, params)
        
        % update transformation model
        obj.Transfo.Params = params;
        
        % transform moving points
        pointsT = transformPoint(obj.Transfo, obj.Points);
        
        % compute metric
        minDist = minDistancePoints(obj.RefPoints, pointsT);
        res = sum(minDist .^ 2);
    end
    
end

end