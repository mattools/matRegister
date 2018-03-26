# matRegister
Matlab library for registration of images and geometric models

MatRegister is a library for Image registration with Matlab. 

The MatRegister library focusses on iconic registration. It follows ITK's registration general workflow, 
and implements the following components:
* several parametric transform models (affine, polynomial, B-Spline...)
* Image-to-image metrics (sum of square differences, mutual information)
* optimizer classes, some of them wrapping the optimization toolbox from matlab.
* various image interpolation strategies (linear, nearest-neighbor)
* image resampling strategies

The MatRegister library rely on the `Image` class, also available on GitHub: https://github.com/mattools/matlab-image-class. 

To install, run the `setupMatRegister.m`, which can be found in the `matRegister` directory.

The library is under development, handle with care...
