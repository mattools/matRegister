classdef SumSquaredDistanceToClosestPoint2D < handle
% SSD metric for two 2D point clouds.
%
%   output = SumSquaredDistanceToClosestPoint2D(REFPTS, PTS, TRANSFO)
%
%   This version computes distances by applying transform to the point set
%   PTS, iterating over transformed points and computing distance to the
%   closest point among reference point set REFPTS.
%
%   Example
%   SumSquaredDistanceToClosestPoint2D
%
%   See also
%     SumSquaredDistanceToClosestPoint3D
 
% ------
% Author: David Legland
% e-mail: david.legland@inrae.fr
% INRAE - BIA Research Unit - BIBS Platform (Nantes)
% Created: 2026-03-24,    using Matlab 9.7.0.1247435 (R2019b) Update 2
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
    function obj = SumSquaredDistanceToClosestPoint2D(refPoints, points, transfo, varargin)
        % Creates a new metric based on point to point distances.
        %
        % Input arguments are the reference point set, the moving set of
        % points, and the transform model.
        %
        % metric = SumSquaredDistanceToClosestPoint2D(refPts, pts, transfo);
        
        if nargin < 3
            error('Require three input arguments');
        end
        
        if  size(refPoints, 2) ~= 2 || size(points, 2) ~= 2
            error('Require 2D points as inputs');
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
    
    function [res, grad] = evaluateWithGradient(obj, params)
        
        % update transformation model
        obj.Transfo.Params = params;
        
        % transform moving points
        pointsT = transformPoint(obj.Transfo, obj.Points);
        
        % identify which reference points are the closest
        % (and keep distance)
        if strcmpi(obj.MatchingMethod, 'Naive')
            inds = findClosestPoint(pointsT, obj.RefPoints);
        elseif strcmpi(obj.MatchingMethod, 'KDTree')
            inds = knnsearch(obj.KDTree, pointsT);
        else
            error('Unknow point matching method');
        end

        refPoints = obj.RefPoints(inds,:);

        diff = pointsT - refPoints;

        % compute metric
        res = sum(sum(diff .^ 2));

        % compute average of gradient of points
        nPoints = size(obj.Points, 1);
        grad = zeros(size(params));
        for i = 1:nPoints
            % retrieve parametric jacobian for moving points
            jac = parametricJacobian(obj.Transfo, obj.Points(i,:));
            % update gradient
            grad = grad + diff(i,:) * jac;
        end
        grad = grad / nPoints;
    end
    
end

end