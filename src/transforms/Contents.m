% TRANSFORMS
%
% A collection of classes for representing paramtric transforms.
%
% Simple (affine) models
%   Translation                    - Defines a translation in ND space
%   Motion2D                       - Composition of a rotation (around origin) and a translation
%   AffineTransform                - Abstract class for AffineTransform
%   TranslationModel               - Transformation model for a translation defined by ND parameters
%   CenteredAffineTransformModel3D - One-line description here, please.
%   CenteredEulerTransform3D       - Transformation model for a centered rotation followed by a translation
%   CenteredMotionTransform2D      - Transformation model for a centered rotation followed by a translation
%   MatrixAffineTransform          - An affine transform defined by its matrix
%
% Basic hierarchy
%   Transform                      - Abstract class for transform
%   ParametricTransform            - Abstract class for parametric transform ND->ND
%   ComposedTransform              - Compose several transforms to create a new transform
%   ComposedTransformModel         - Compose several transforms, the last one being parametric 
%   CenteredTransformAbstract      - Add center management to a transform
%   InitializedTransformModel      - InitializedTransformModel Encapsulate a parametric and an initial transform
%
% Spline models
%   BSplines                       - Contains several static functions for manipulation of cubic splines
%   CenteredQuadTransformModel3D   - CenteredQuadTransformModel3D  One-line description here, please.
%
