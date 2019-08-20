function drawTransformedGrid(transfo, lx, ly, varargin)
% Draw the result of a transform applied to a grid.
%
%   usage:
%   drawTransformedGrid(T, LX, LY)
%   T: an instance of Transform class
%   LX, LY: specify spacing of grid vertices, see meshgrid for details.
%
%   drawTransformedGrid(T, LX, LY, NPTS)
%   Speccify the number of line segments used to represent the junction
%   between two grid vertices. Default is 5.
%
%   drawTransformedGrid(..., PNAME, PVALUE)
%   Specify additional parameters to the plot function, using parameter
%   name-value pairs.
%
%   Example
%     % Creates a basic motion transform, and disply the result of grid
%     % stransformation over a default image.
%     T = CenteredMotionTransform2D([10 20 30]);
%     setCenter(T, [128 128]);
%     img = imread('cameraman.tif');
%     figure; imshow(img); hold on;
%     drawTransformedGrid(T, 16:32:256, 16:32:256, 'color', 'g')
%
%   See also
%     Transform
 
% ------
% Author: David Legland
% e-mail: david.legland@inra.fr
% Created: 2018-08-09,    using Matlab 9.4.0.813654 (R2018a)
% Copyright 2018 INRA - Cepia Software Platform.

nSteps = 5;
if ~isempty(varargin) && isscalar(varargin{1})
    nSteps = varargin{1};
    varargin(1) = [];
end

if length(varargin) ~= 1
    varargin = [{'color', 'g'} varargin];
end

% compute sub-sampled bases
lx2 = lx(1):((lx(2)-lx(1))/nSteps):lx(end);
ly2 = ly(1):((ly(2)-ly(1))/nSteps):ly(end);

% x-grid
[x, y] = meshgrid(lx2, ly);
pts = [x(:) y(:)];
pts2 = transformPoint(transfo, pts);
x2 = reshape(pts2(:,1), size(x))';
y2 = reshape(pts2(:,2), size(x))';
plot(x2, y2, varargin{:});

hold on;

% y-grid
[x, y] = meshgrid(lx, ly2);
pts = [x(:) y(:)];
pts2 = transformPoint(transfo, pts);
x2 = reshape(pts2(:,1), size(x));
y2 = reshape(pts2(:,2), size(x));
plot(x2, y2, varargin{:});

