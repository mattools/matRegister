%DEMO_BSPLINETRANSFORMMODEL2D_CAMERAMAN  One-line description here, please.
%
%   output = demo_BSplineTransformModel2D_cameraman(input)
%
%   Example
%   demo_BSplineTransformModel2D_cameraman
%
%   See also
%
 
% ------
% Author: David Legland
% e-mail: david.legland@inra.fr
% Created: 2018-08-09,    using Matlab 9.4.0.813654 (R2018a)
% Copyright 2018 INRA - Cepia Software Platform.


%% Read image and create transform

img = imread('cameraman.tif');

transfo = BSplineTransformModel2D([3 3], [64 64], [64 64]);
% init parameters to arbitrary perturbation
transfo.Params = [...
    +1 +1   -1 +1    0  0 ...
    +1 -1   -1 -1   +1 -1 ...
     0  0   -1 +1   +1 +1 ...
    ] * 30;

figure(1); close;
figure(1); imshow(img);
hold on; drawGrid(transfo);
drawVertexShifts(transfo, 'm');


%% points within the whole image

lx = 10:10:250;
ly = 10:10:250;
[x, y] = meshgrid(lx, ly);
pts = [x(:) y(:)];

pts2 = transformPoint(transfo, pts);
drawPoint(pts2, 'g.');

drawTransformedGrid(transfo, 10:10:250, 10:10:250)


%% Transform image

% create resampling grid
[x, y] = meshgrid(1:256, 1:256);
pts = [x(:) y(:)];

% compute transformed positions
pts2 = transformPoint(transfo, pts);

% create new image
img2 = reshape(imEvaluate(img, pts2), size(img));

% display result
figure;
imshow(img2, [0 255]);


%% Compute Map of determinant of Jacobian

% compute the 2-by-2-by-Npts jacobian matrix
jac = jacobianMatrix(transfo, pts);

% compute map of jacobian determinant
jacMap = zeros(size(img));
for i = 1:length(pts)
    jacMap(i) = det(jac(:,:,i));
end

% use log scale for symmetric behaviour
logJacMap = log2(jacMap);

% display with color coding
figure;
imshow(logJacMap, [-2 2]); 
colormap (blue2White2Red); 

% add color bar, using log values
hcb = colorbar;
ticks = get(hcb, 'Ticks');
labels = strtrim(cellstr(num2str(power(2, ticks)', '%.2g')));
set(hcb, 'TickLabels', labels)

