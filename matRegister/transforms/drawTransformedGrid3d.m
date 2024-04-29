function drawTransformedGrid3d(transfo, lx, ly, lz, varargin)
% Draw the result of a transform applied to a grid.
%
%   Usage:
%   drawTransformedGrid3d(T, LX, LY, LZ)
%   T: an instance of Transform class
%   LX, LY: specify spacing of grid vertices, see meshgrid for details.
%
%   drawTransformedGrid3d(T, LX, LY, LZ, NPTS)
%   Specify the number of line segments used to represent the junction
%   between two grid vertices. Default is 5.
%
%   drawTransformedGrid3d(..., PNAME, PVALUE)
%   Specify additional parameters to the plot function, using parameter
%   name-value pairs.
%
%   Example
%     % Creates a basic motion transform, and display the result of grid
%     % transformation over a default image.
%     T = CenteredMotionTransform2D([10 20 30]);
%     setCenter(T, [128 128]);
%     img = imread('cameraman.tif');
%     figure; imshow(img); hold on;
%     drawTransformedGrid3d(T, 16:32:256, 16:32:256, 'color', [.5 .5 .5])
%
%   See also
%     Transform, meshgrid, drawTransformedGrid
%
 
% ------
% Author: David Legland
% e-mail: david.legland@inrae.fr
% Created: 2018-08-09,    using Matlab 9.4.0.813654 (R2018a)
% Copyright 2018 INRA - Cepia Software Platform.

nSteps = 5;
if ~isempty(varargin) && isscalar(varargin{1})
    nSteps = varargin{1};
    varargin(1) = [];
end

if length(varargin) ~= 1
    varargin = [{'color', 'b'} varargin];
end

% compute sub-sampled bases
lx2 = lx(1):((lx(2)-lx(1))/nSteps):lx(end);
ly2 = ly(1):((ly(2)-ly(1))/nSteps):ly(end);
lz2 = lz(1):((lz(2)-lz(1))/nSteps):lz(end);

% x-grid
[x, y, z] = meshgrid(lx2, ly, lz);
pts = [x(:) y(:) z(:)];
pts2 = transformPoint(transfo, pts);
dims2 = [length(lx2), length(ly) * length(lz)];
x2 = reshape(permute(reshape(pts2(:,1), size(x)), [2 1 3]), dims2);
y2 = reshape(permute(reshape(pts2(:,2), size(x)), [2 1 3]), dims2);
z2 = reshape(permute(reshape(pts2(:,3), size(x)), [2 1 3]), dims2);
plot3(x2, y2, z2, varargin{:});

hold on;

% y-grid
[x, y, z] = meshgrid(lx, ly2, lz);
pts = [x(:) y(:) z(:)];
pts2 = transformPoint(transfo, pts);
dims2 = [length(ly2), length(lx) * length(lz)];
x2 = reshape(permute(reshape(pts2(:,1), size(x)), [1 2 3]), dims2);
y2 = reshape(permute(reshape(pts2(:,2), size(x)), [1 2 3]), dims2);
z2 = reshape(permute(reshape(pts2(:,3), size(x)), [1 2 3]), dims2);
plot3(x2, y2, z2, varargin{:});


% z-grid
[x, y, z] = meshgrid(lx, ly, lz2);
pts = [x(:) y(:) z(:)];
pts2 = transformPoint(transfo, pts);
dims2 = [length(lz2), length(lx) * length(ly)];
x2 = reshape(permute(reshape(pts2(:,1), size(x)), [3 1 2]), dims2);
y2 = reshape(permute(reshape(pts2(:,2), size(x)), [3 1 2]), dims2);
z2 = reshape(permute(reshape(pts2(:,3), size(x)), [3 1 2]), dims2);
plot3(x2, y2, z2, varargin{:});

