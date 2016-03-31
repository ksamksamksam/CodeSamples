/****************************************************************
 * @file
 * @author Samantha Raja <ksam@seas.upenn.edu>									
 * @date Summer 2010
 *
 * @section DESCRIPTION 
 * This file contains the struct definitions for light and the class declaration of the 
 * RayTracer class, which contains scene information and renders a raytraced bitmap image.
 *
 ****************************************************************/

#ifndef RAYTRACER_H
#define RAYTRACER_H

#include "object.h"

#define SMALL .00000001 //close enough to zero to treat as zero

/****************************************************************/

//This struct defines the light properties
struct Light {
	vec3 pos;
	vec3 color; //color of the light
};

/****************************************************************/

class RayTracer{
public:
	//Constructors and destructor
	RayTracer();
	RayTracer(int w, int h, int md);
	RayTracer(const RayTracer& r);
	~RayTracer();

	//Public methods
	void addObject(int type, vec3 pos, double scale, Material mat, vec3 axis, double angle);
	void addLight(vec3 pos, vec3 color);

	void setDefaultScene();

	void writeImage(); // "go" button, basically

	//Public member variables
	int height;
	int width;

protected:
	//Protected methods
	vec3 traceRay(vec3 start, vec3 ray, int depth);

	double planeCollision(vec3 origin, vec3 direction, Object* p);
	double sphereCollision(vec3 origin, vec3 direction, Object* s);
	double intersectTest(vec3 start, vec3 ray, Object* obj);
	vec3 refract(vec3 ray, vec3 normal, double IOR1, double IOR2);
	vec3 clamp(vec3 v);

	//Member variables
	std::vector<Object*> objects; //list of objects in the scene
	std::vector<Light*> lights;
	int maxDepth; //the maximum depth of recursion of the ray tree
};

#endif