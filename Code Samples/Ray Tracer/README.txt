NAME: Samantha Raja
INTERVIEWING FOR: Technical Artist


CODE SAMPLE


DESCRIPTION

This project is a simple raytracer that I wrote. The pipeline is as follows:  main method->set up scene->build image->trace a ray through each pixel->write image to file

The application processes two primitives: spheres and planes.  Each object contains a Material which specifies diffuse, ambient, specular, refractive, and reflective weights.  It also specifies a color.  Light calculations are made based on the Phong local illumination model, with an additional term for reflective color.  The application can also process multiple lights of arbitrary colors.

The objects and materials in the scene can be changed in the setDefaultScene() method in the RayTracer class.  Currently it is a Cornell box containing two spheres.

The code was written in 2010 in C++.


CLASS AND STRUCT DESCRIPTIONS

Material: This struct contains a set of shading properties.

Light: This struct contains a set of light properties.

Object: An object can be of two primitive types--sphere or plane--and contains a Material object, a scale, an axis and angle or rotation, and a position.

RayTracer: The muscle.  This class contains the traceRay method and all of its helper methods, including intersection tests, a scene setup, and a writeImage method.


EXTERNAL LIBRARIES

EasyBMP: This is an image processing library that I used to write the pixel information to an image file.  It was created by Paul Macklin of the EasyBMP Project. <http://easybmp.sourceforge.net/>

algebra3.h: This is a version of the vector and matrix math library found in NVidia's "Graphics Gems IV" which has been modified by J. Nagle.

Both libraries are open source.