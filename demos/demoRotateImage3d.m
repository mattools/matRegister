function varargout = demoRotateImage3d(varargin)
%DEMOROTATEIMAGE3D  One-line description here, please.
%
%   output = demoRotateImage3d(input)
%
%   Example
%   demoRotateImage3d
%
%   See also
%
%
% ------
% Author: David Legland
% e-mail: david.legland@grignon.inra.fr
% Created: 2010-07-27,    using Matlab 7.9.0.529 (R2009b)
% Copyright 2010 INRA - Cepia Software Platform.


%% Lecture

% lecture de l'image de test
imgDir = fullfile('d:', 'images', 'amib', 'clermont', 'roi_memri');
fileName = fullfile(imgDir, 'BRNOR39e5p1', 'BRNOR39e5p1Ball5Msk160.mhd');
img = metaImageRead(fileName);

% selectionne la zone interessante, et convertit en 256 niveaux de gris
imgCrop = img(1:160, 31:130, 81:160);
imgCrop = imRescale(imgCrop);

figure(1); clf; hold on;
orthoSlices(imgCrop);


%% Creation du maillage

% binarisation et nettoyage
imgCropBin = imopen(imgCrop>0, ones([3 3 3]));

% creation de la transformation a appliquer
rot = createEulerAnglesRotation(10*pi/180, 20*pi/180, 30*pi/180);
rot2 = recenterTransform3d(rot, [50 80 40]);

% calcul de la surface du cerveau
[f0 v0] = isosurface(imgCropBin);

% affiche la surface
figure(2); clf; hold on;
drawMesh(v0, f0, 'facecolor', 'g', 'linestyle', 'none');

% calibre affichage
set(gcf, 'renderer', 'opengl')
camlight left
camlight right
lighting phong

% transforme les sommets du maillage
vt = transformPoint3d(v0, rot2);

% affiche le maillage transforme
drawMesh(vt, f0, 'facecolor', 'b', 'linestyle', 'none');


%% Interpolation de l'image

I = Image3D(imgCrop);

% definition de l'interpolateur
interp = LinearInterpolator3D(I);

samp = ImageResampler(I);

% Definit la fonction de transfo inverse de l'image
trans = MatrixAffineTransform(inv(rot2));
tim = BackwardTransformedImage(interp, trans);

% echantillonnage de l'image transformee
I2 = samp.resample(tim);

% binarise l'image transformee pour comparer les contours
imgTransBin = imopen(I2.getBuffer>0, ones([3 3 3]));

% calcul isosurface du resultat
[ft vt] = isosurface(imgTransBin, .5);

% affiche en surimpression, les maillages violet et bleu doivent se
% superposer
drawMesh(vt, ft, 'facecolor', 'm', 'linestyle', 'none');

