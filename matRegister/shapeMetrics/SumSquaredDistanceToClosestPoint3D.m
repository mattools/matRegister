classdef SumSquaredDistanceToClosestPoint3D < handle
% SSD metric for two 3D point clouds.
%
%   output = SumSqDistancesVertexToClosestVertex3d(PTS1, PTS2, TRANSFO)
%
%   This version computes distances by applying transform to the point set
%   2, iterating over transformed points and computing distance to closest
%   point among reference point set 1.
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
    
    % The name of the matching method
    % Must be one of "naive", "kdtree"
    MatchingMethod = 'KDTree';

    % If matching method is "kd-tree", stores the KDTree for repeated queries.
    KDTree;
end

%% Constructor method
methods
    function obj = SumSquaredDistanceToClosestPoint3D(refPoints, points, transfo, varargin)
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
        
        % Process optional arguments
        while length(varargin) > 1
            pname = varargin{1};
            if ~ischar(pname)
                error('Requires paramaters given as name-value pairs');
            end
            if strcmpi(pname, 'matchingMethod')
                obj.MatchingMethod = varargin{2};
            else
                error("Unknown option name: %s", pname);
            end

            varargin(1:2) = [];
        end


        if strcmpi(obj.MatchingMethod, 'kdTree')
            obj.KDTree = KDTreeSearcher(obj.RefPoints);
        end

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
        if strcmpi(obj.MatchingMethod, 'Naive')
            minDist = minDistancePoints(pointsT, obj.RefPoints);
        elseif strcmpi(obj.MatchingMethod, 'KDTree')
            [~, minDist] = knnsearch(obj.KDTree, pointsT);
        else
            error('Unknow point matching method');
        end

        res = sum(minDist .^ 2);
    end
    
end

end