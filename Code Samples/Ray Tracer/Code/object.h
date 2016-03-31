/****************************************************************
 * @file                                                            
 * @author  Samantha Raja <ksam@seas.upenn.edu>                                 
 * @date  Summer 2010
 *
 * @section DESCRIPTION
 * This file contains the Object class, which represents the two 
 * types of primitives my application processes: spheres and planes.  
 * It also contains the definition for the Material class, 
 * which represents a set of shading properties.
 *
 ****************************************************************/

#ifndef OBJECT_H
#define OBJECT_H

#include <vector>
#include "algebra3.h" //vector and matrix math library

/****************************************************************/
//MATERIAL
//Contains a set of shading properties.
struct Material {
public:
    Material( double Kd, double Ks, int Kn, double Ka, double Kt, double Kr, double Kref, vec3 color ) 
        : Kd(Kd), Ks(Ks), Kn(Kn), Ka(Ka), Kt(Kt), Kr(Kr), Kref(Kref), color(color) {}

    double Kd; //diffuse weight
    double Ks; //specular weight
    int Kn; //specular Phong exponent
    double Ka; //ambience weight
    double Kt; //transmittance
    double Kr; //index of refraction
    double Kref; //mirror reflectiveness of surface
    vec3 color; //color, represented as rgb values from 0-1
};

/****************************************************************/
//OBJECT
//This is the primitive class that defines the objects that can populate a raytraced scene.  
//Each Object instance can be of type SPHERE, CUBE, or PLANE. 

enum {SPHERE, PLANE};
class Object;

//Declaration
class Object{
public:
    //Constructors and destructor
    Object();
    Object(int t, vec3 p, double s, Material* m, vec3 a, double theta);
    Object(const Object& o); 
    ~Object();
    
    //Member variables
    int type; //type of object (i.e. Sphere or Plane)
    vec3 pos; //position of object
    double scale; //scale of the object. 
                  //If SPHERE, this is the radius.  If PLANE, this is the length of one side.
    Material* mat; //material of object
    vec3 axis; //axis of rotation of object, irrelevant if SPHERE
    double angle; //angle of rotation of object, irrelevant if SPHERE
};


//Implementation

//Constructors
inline Object::Object()
{
}

inline Object::Object(const Object& o) : type(type), pos(pos), scale(scale), axis(axis), angle(angle)
{
    mat = new Material(*(o.mat));
}

inline Object::Object(int t, vec3 p, double s, Material* m, vec3 a, double theta)
{
    type = t;
    pos = p;
    scale = s;
    mat = m;
    axis = a;
    angle = theta;
}

//Destructor
inline Object::~Object()
{
    delete mat;
}

/****************************************************************/

#endif OBJECT_H