% TRANSFORMS
%
% A collection of classes for representing parametric transforms.
%
% Non-Parametric Affine transforms (for validation purpose)
%   Transform                      - Abstract class for transform
%   Translation                    - Defines a translation in ND space
%   Motion2D                       - Composition of a rotation (around origin) and a translation
%   AffineTransform                - Abstract class for AffineTransform
%   MatrixAffineTransform          - An affine transform defined by its matrix
%   ComposedTransform              - Compose several transforms to create a new transform
%
% Affine Transform models
%   ParametricTransform            - Abstract class for parametric transform ND->ND
%   MotionModel2D                  - Transformation model for a centered rotation followed by a translation
%   SimilarityModel2D              - Transformation model for a similarity: rotation+scaling followed by a translation
%   AffineTransformModel2D         - Transformation model for a 2D affine transform
%   TranslationModel               - Transformation model for a translation defined by ND parameters
%   CenteredAffineTransformModel3D - Transformation model for a centered 3D affine transform followed by a translation
%   CenteredEulerTransform3D       - Transformation model for a centered 3D rotation followed by a translation
%   CenteredMotionTransform2D      - Transformation model for a centered rotation followed by a translation
%
% Polynomial and spline models
%   BSplineTransformModel2D        - Cubic B-Spline Transform model in 2D
%   QuadPolynomialTransformModel2D - One-line description here, please.
%   CenteredQuadTransformModel3D   - Polynomial 3D transform up to degree 2 (3*10=30 parameters)
%   CenteredQuadTransformModel2D   - Polynomial 2D transform up to degree 2 (2*6=12 parameters)
%   BSplines                       - Contains several static functions for manipulation of cubic splines
%   RadialScalingTransform2D       - Radial scaling transform in 2D
%
% Utility classes for parametric transforms
%   ComposedTransformModel         - Compose several transforms, the last one being parametric 
%   CenteredTransformAbstract      - Add center management to a transform
%   InitializedTransformModel      - Encapsulation of a parametric and an initial transform
%
% Utility methods
%   drawTransformedGrid            - Draw the result of a transform applied to a grid

