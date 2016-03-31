/****************************************************************
 * @file                                                            
 * @author  Samantha Raja <ksam@seas.upenn.edu>                                 
 * @date  Summer 2010
 *
 * @section DESCRIPTION
 * This file contains the implementation of the RayTracer class, 
 * which contains scene information and renders a raytraced image.
 *
 ****************************************************************/

#include "easybmp.cpp" //BMP image writing library
#include "raytracer.h"

#define air_IOR 1.0003 //index of refraction of air is 1.0003
#define OUTPUT_IMAGE "spec_mirror_02.jpg" //name of the image file we want to write to

BMP image; //output image map

/****************************************************************/
//CONSTRUCTORS

RayTracer::RayTracer()
{
    image.SetSize(512,512);
    maxDepth = 3;
}

RayTracer::RayTracer(int w, int h, int md)
{
    width = w;
    height = h;
    image.SetSize(width, height);
    maxDepth = md;
}

RayTracer::RayTracer(const RayTracer& r)
{
    image.SetSize(r.width, r.height);
    maxDepth = r.maxDepth;

    for (unsigned int i = 0; i < r.lights.size(); i++)
        lights.push_back(new Light(*(r.lights[i])));

    for (unsigned int i = 0; i < r.objects.size(); i++)
        objects.push_back(new Object(*(r.objects[i])));
}

/***************************************************************/
//PUBLIC METHODS

/**
 * addObject
 *
 * Creates an object with the given properties and adds it to the object list.
 *
 * @param  type Type of object (either SPHERE or PLANE)
 * @param  pos Position of the object
 * @param  scale Scale of the object
 * @param  mat Material belonging to the object
 * @param  axis Axis of rotation
 * @param  angle Angle of rotation
 * 
 */
void RayTracer::addObject(int type, vec3 pos, double scale, Material mat, vec3 axis, double angle)
{
    Object* obj = new Object();
    obj->scale = scale;
    obj->type = type;
    obj->pos = pos;
    obj->mat = new Material(mat);
    obj->axis = axis;
    obj->angle = angle;

    objects.push_back(obj);
}

/**
 * addLight
 *
 * Add a light with the given properties to the scene.
 *
 * @param  pos Position of the light
 * @param  color Color of the light
 *
 */
void RayTracer::addLight(vec3 pos, vec3 color)
{
    Light* l = new Light();
    l->pos = pos;
    l->color = color;
    lights.push_back(l);
}

/**
 * setDefaultScene
 *
 * Sets up a simple Cornell Box scene with one red wall, one green wall, a blue-green specular sphere, and a mirror sphere.
 * This is where one would go to manipulate the scene and try different combinations.
 *
 */
void RayTracer::setDefaultScene()
{
    //Define our materials
    Material redDiffuse(1,0,0,.1,0,0,0, vec3(1,0,0)); //material for red wall
    Material greenDiffuse(1,0,0,.1,0,0,0,vec3(0,1,0)); //material for green wall
    Material whiteDiffuse(1,0,0,.1,0,0,0,vec3(1,1,1)); //material for white wals
    Material mirrorReflect(0,0,0,0,0,0,1,vec3(1,0,0)); //mirror material
    Material glass(0,.5,0,0,1,1.52,.1,vec3(1,1,1)); //glass material
    Material blueSpec(1,.7,12,.1,0,0,.12,vec3(.2,.6,1)); //shiny blue material
    Material yellowDiffuse(1,0,0,.1,0,0,0, vec3(1,1,0)); //yellow diffuse

    //Add the spheres
    addObject(SPHERE, vec3(120,150,-750), 9000, mirrorReflect, vec3(0,1,0), 0); //glass sphere
    addObject(SPHERE, vec3(-120,-100,-850), 9000, blueSpec, vec3(0,1,0), 0); //blue shiny sphere

    //Make the four walls of the box
    addObject(PLANE, vec3(-250,0,-750), 500, redDiffuse, vec3(0,1,0), 90); //left wall
    addObject(PLANE, vec3(250,0,-750), 500, greenDiffuse, vec3(0,1,0), -90); //right wall
    addObject(PLANE, vec3(0,-250,-750), 500, whiteDiffuse, vec3(1,0,0), -90); //bottom wall
    addObject(PLANE, vec3(0,250,-750), 500, whiteDiffuse, vec3(1,0,0), 90); //top wall
    addObject(PLANE, vec3(0,0,-999), 500, whiteDiffuse, vec3(0,1,0), 0); //back wall

    //Add lights
    addLight(vec3(0,-200,-750), vec3(1,1,1));
//  addLight(vec3(0,200,-750), vec3(.5,.5,.5)); //second light
}

/**
 * writeImage
 *
 * Writes the raytraced image to the file defined by OUTPUT_IMAGE and prints out its progress as it goes along.
 *
 */
void RayTracer::writeImage()
{
    //record progress
    double percent = 0;

    //for each screen pixel, shoot a ray from the centered camera through the pixel and trace the ray.
    for (int x = -250; x < 250; x++)
    {
        percent += .2;

        //every 2 percent or so, print the total percentage done
        if (int(percent*10)%10 == 0)
        {
            cout<<percent<<" percent done.\n";
        }
        for (int y = -250; y < 250; y++)
        {
            vec3 pix = vec3(x,y,-500);

            //trace ray through pixel
            vec3 color = traceRay(vec3(0,0,0), pix, 0);

            //record resulting color to the image
            image(x+250,y+250)->Red = color[0]*255;
            image(x+250,y+250)->Green = color[1]*255;
            image(x+250,y+250)->Blue = color[2]*255;
        }
    }

    //write the image to a file
    image.WriteToFile(OUTPUT_IMAGE);
}

/**
 * traceRay
 *
 * This is the meat of our program.  Takes a start position, a direction vector, and a depth number as input
 * and outputs the color value at the final resting point of the ray.  The method is only allowed to recurse maxDepth times,
 * at which point it returns black as the output color.
 *
 * @param  start Origin of our ray
 * @param  ray Direction vector of our ray
 * @param  depth Number of times we have recursively called traceRay
 *
 */
vec3 RayTracer::traceRay(vec3 start, vec3 ray, int depth)
{
    vec3 reflectedRay, refractedRay, transmittedRay; //ray vectors we will use
    vec3 reflectedColor, refractedColor;    //color vectors we will pass into recursive calls to traceRay
    vec3 refl, refr, outColor;  //final color vectors we will use in the shading equation

    //if we have already recursively called traceRay more than maxDepth times, stop and set outColor to black
    if (depth > maxDepth)
    {
        outColor = vec3(0,0,0);
        return outColor;
    }

    ray.normalize(); //make sure our ray is unit length

    //CHECK FOR INTERSECTIONS
    vec3 intersectPt;
    double minResult = -1; //Minimum distance to an intersected object.  Set to -1 initially.
    Object* intersectObj; //We want this to point to the closest intersected object.

    //Traverse through our object list and check for intersections.  If an object is intersected, check to see if the distance
    //from our start point to the object is less than the smallest intersection distance so far.  If so, this distance
    //shall be the new smallest distance.  The final minResult value should reflect the minimum distance and intersectObj
    //should point to the closest object.
    for (unsigned int i = 0; i < objects.size(); i++)
    {
        Object* obj = objects[i];
        double result = intersectTest(start, ray, obj);

        if (result > 0)
        {
            //if our result is less than minResult, or if minResult hasn't been set yet, set minResult to result and set
            //intersectObj to obj.
            if (result - minResult < 0 || minResult == -1)
            {
                minResult = result;
                intersectObj = obj;
            }
        }
    }

    //If the ray does not intersect any objects, set outColor to black and exit the method
    if (minResult == -1)
    {
        outColor = vec3(0,0,0);
        return outColor;
    }

    intersectPt = start + minResult * ray; //calculate point of intersection

    vec3 normal;

    if (intersectObj->type != PLANE)
        normal = (intersectPt-intersectObj->pos).normalize(); //calculate the normal.  
                                                              //normal = point of intersection - center of object

    //Find the normal of a plane by rotating a viewer-facing plane normal by the axis-angle rotation
    else
        normal = (rotation3D(intersectObj->axis, intersectObj->angle) * vec3(0,0,1)).normalize();
    
    Material* objMat = intersectObj->mat;
    
    //CALCULATE REFLECTED COLOR
    if (objMat->Kref > 0) //if the object's material has specular reflectivity
    {
        //Rr = Ri - 2*N*(Ri*N)
        reflectedRay = ray - 2 * normal * ray * normal;

        //call traceRay on the reflected direction with the intersection point as the origin 
        //and store that color value in reflectedColor
        reflectedColor = traceRay(intersectPt, reflectedRay, depth+1);
        refl = objMat->Kref * reflectedColor; //refl is the weighted reflected color
    }

    else 
        refl = vec3(0,0,0); //else, set the reflected color to black

    //CALCULATE REFRACTED COLOR
    if (objMat->Kt > 0) //if the object's material has transmittance
    {
        double obj_IOR = objMat->Kr;
        double refrRatio = air_IOR / obj_IOR; //take the ratio of the refractive index of air to the object material
        
        //First we calculate the refracted ray INTO the surface
        vec3 refrIn = refract(ray, normal, air_IOR, obj_IOR);

        //if the ray is at the critical angle, set the refracted color to black
        if (refrIn == NULL)
        {
            refr = vec3(0,0,0);
        }

        else{
            refrIn.normalize();

            //Now, we want to find the exit point of the ray--where it intersects with the object from the inside.  But, our
            //intersection tests ignore rays that pass through the inside of the object, so we need to take a different approach.
            //We will shoot a ray in the opposite direction from outside the object.  The point of intersection will be the same
            //as if we had shot the ray from the inside.

            vec3 outsidePt = intersectPt + 2 * intersectObj->scale * refrIn;  //This is our outside point.  We know it must be outside
                                                                          //the object because its distance from the intersection 
                                                                          //point is twice the object's scale

            double refrDist = intersectTest(outsidePt, -refrIn, intersectObj);

            if (refrDist == -1) //If the object is flat and there is no second intersection, use the ray into the surface
            {
                refractedRay = refrIn;
            }

            else{
                vec3 secondIntersectPt = outsidePt + refrDist * (-refrIn); //compute second intersect point

                vec3 normal2 = (intersectObj->pos - secondIntersectPt).normalize(); //compute normal from the intersect point to
                                                                                    //the object's center

                refractedRay = refract(refrIn, normal2, obj_IOR, air_IOR); //Now we can find the refracted ray OUT of the object
            }

            //if our first refracted ray was an critical angle, set the refracted color to black
            if (refractedRay == NULL)
            {
                refr = vec3(0,0,0);
            }

            else{
                //call traceRay on the refracted direction with the intersection point as the origin 
                //and store that color value in refractedColor
                refractedColor = traceRay(intersectPt, refractedRay, depth+1);
                refr = objMat->Kt * refractedColor; //refr is the weighted refracted color
            }
        }
    }

    else 
        refr = vec3(0,0,0); //else, set the refr color to black

    //CALCULATE OVERALL COLOR

    //first add the ambient and refracted colors--these are unaffected by shadows or light source
    outColor = objMat->Ka * objMat->color + refr;

    //for each light, calculate the diffuse and specular terms and sum all of them together
    for (unsigned int i = 0; i < lights.size(); i++)
    {
        vec3 lightPos = lights[i]->pos;
        vec3 lightNorm = (lightPos-intersectPt).normalize(); //normalized vector from intersection point to light source
        bool shadowCheck = false; //we store the result of our shadow check here

        //shadow check
        for (unsigned int j = 0; j < objects.size(); j++)
        {
            //shoot rays toward light source and look for intersections
            double t = intersectTest(intersectPt, lightNorm, objects[j]);

            //if there is an intersection and it is not a self-intersection or an intersection with an object
            //that lies beyond the light, set shadowCheck to true
            if (t > SMALL && (t * lightNorm).length() < (lightPos-intersectPt).length())
            {
                shadowCheck = true;
            }
        }

        //If there are no objects blocking the path from the intersect point to the light source, the diffuse and specular
        //contribution from the light source
        if (shadowCheck == false)
        {
            //calculate diffuse color 

            //light attenuation factor: light falls off at a rate proportional to the inverse square of the distance 
            //between surface and source light [Lambert's first law].  1/(c1+c2*d+c3*d^2).  Equation proposed by Foley and Van Dam.
            double col1 = 1/(1 + .0007 * ((lightPos-intersectPt).length())+ .000007 * ((lightPos-intersectPt).length() *
                (lightPos - intersectPt).length())); 
            
            double cos = (lightNorm) * normal; //brightness depends on angle between light source direction and surface 
                                               //normal [Lambert's second law]. cos_theta = L * N.
            
            //clamp negative angle values to 0
            if (cos < 0) 
                cos = 0;

            double fac = cos * col1; //combine into one brightness factor

            //Diffuse term = color * Kd * brightness_factor
            vec3 diffColor = objMat->color*(objMat->Kd * fac);
            diffColor = clamp(diffColor); //clamp color values to [0,1]

            //calculate specular color

            //Specular term = Ks*((light vector reflected about normal * viewing vector)^Kn)
            vec3 specColor = fac * (objMat->Ks * pow((2 * (normal * -ray) * (normal * lightNorm) - 
                (lightNorm * -ray)), objMat->Kn));
            specColor = clamp(specColor); //clamp color values to [0,1]

            //add them to the overall color and reflection color and multiply by the light color
            outColor += vec3(diffColor[0] * lights[i]->color[0], diffColor[1] * lights[i]->color[1], 
                diffColor[2] * lights[i]->color[2])  + vec3(specColor[0] * lights[i]->color[0], 
                specColor[1] * lights[i]->color[1], specColor[2] * lights[i]->color[2]) + refl;
            outColor = clamp(outColor);  //clamp color values to [0,1]
        }
    }

    //clamp color values to [0,1]
    outColor = clamp(outColor);
    return outColor;
}


/***********************************************************************************/
//PROTECTED METHODS

/** 
 * sphereCollision
 *
 * Calculates the intersection point of a ray with the given origin and direction
 * with the given Sphere.
 *
 * @param  origin Origin of ray
 * @param  direction Direction of ray
 * @param  s Collision sphere
 * 
 * @return  double p such that [ origin + p * direction.normalized = point of intersection ] if an intersection exists
 *          -1 if there is no intersection
 *
 */
double RayTracer::sphereCollision(vec3 origin, vec3 direction, Object* s)
{
    direction.normalize(); //make sure the direction vector is unit sized
    vec3 dst = origin - s->pos; //distance vector between origin of ray and sphere

    double B = dst * direction; 
    double C = dst * dst - s->scale;
    double D = B * B - C;
    double p = D > 0 ? -B - sqrt(D) : -1; //distance along ray that intersection point lies if it exists; -1 if it does not

    vec3 n((origin + p * direction) - s->pos); //normal out of intersection point

    // If the normal and the ray are orthogonal or the origin of the ray 
    // lies within the sphere, not considered an intersection
    if (n * direction == 0 || (origin - s->pos).length()- sqrt(s->scale) < SMALL)
    {
        return -1;
    }

    return p;
}

/** 
 * planeCollision
 *
 * Calculates the intersection point of a ray with the given origin and direction
 * with the given Plane.
 *
 * @param  origin Origin of ray
 * @param  direction Direction of ray
 * @param  p Collision plane
 * 
 * @return  double p such that [ origin + p * direction.normalized = point of intersection ] if an intersection exists
 *          -1 if there is no intersection
 *
 */
double RayTracer::planeCollision(vec3 origin, vec3 direction, Object* p)
{
    //the rotation matrix of the plane
    mat4 planeRot = rotation3D(p->axis, p->angle);

    //pick an arbitrary point on the plane.  We found a point that lies on an edge of the plane, then rotated it to reflect
    //the plane's rotation.
    vec3 arbPt = p->pos + planeRot * vec3(p->scale / 2.0,0,0);

    //Find the normal by rotating a viewer-facing plane normal by the appropriate amount
    vec3 norm = (planeRot * vec3(0,0,1)).normalize();

    //Check if ray is parallel to or in the plane.  If so, return -1
    if (direction * norm == 0)
        return -1;

    //lineDist is the distance along the ray at which the intersection occurs.
    double lineDist = ((arbPt-origin) * norm)/(direction * norm);
    
    //check to make sure the intersection point lies within the plane.  For simplicity, we will do this by rotating the intersection point
    //into the plane's local coordinates and making sure it lies within the boundaries.

    vec3 intersectPt = ((origin+lineDist * direction) - p->pos); //intersection point
    
    //rotation matrix
    mat4 rayRot = rotation3D(p->axis, -p->angle);

    //rotated point of intersection
    vec3 rotatedRay = rayRot * intersectPt;

    //if the point lies outside any of the 4 lines that make up the plane's edges, return -1
    double posY = p->pos[1];
    double posX = p->pos[0];
    double boundCheck = p->scale;
    if (rotatedRay[0] > posX + boundCheck || rotatedRay[0] < posX - boundCheck || rotatedRay[1] > posY + boundCheck 
        || rotatedRay[1] < posY - boundCheck)
        return -1;

    return lineDist;
}

/** 
 * intersectTest
 *
 * Checks whether the given ray intersects with the given object.  If it does, return the intersection distance; if not
 * return -1.
 *
 * @param  start Origin of ray
 * @param  ray Direction of ray
 * @param  obj We are checking this object for an intersection
 * 
 * @return  double p such that [ origin + p * direction.normalized = point of intersection ] if an intersection exists
 *          -1 if there is no intersection
 *
 */
double RayTracer::intersectTest(vec3 start, vec3 ray, Object* obj)
{
    double result = -1;

    //Check the object's type and perform the appropriate intersection test.
    if (obj->type == SPHERE)
        result = sphereCollision(start, ray, obj);
    else
        result = planeCollision(start, ray, obj);

    if (result < SMALL)
        result = -1;

    return result;
}

/** 
 * refract
 *
 * Calculates the refracted ray of the given ray at the given normal, through substances having the given indices of refraction,
 *
 * @param  ray Ray that we want to refract
 * @param  normal Normal which we are using as our refraction pivot
 * @param  IOR1 The index of refraction of the first substance
 * @param  IOR2 The index of refraction of the second substance
 * 
 * @return  The refracted ray
 *
 */
vec3 RayTracer::refract(vec3 ray, vec3 normal, double IOR1, double IOR2)
{
    double refrRatio = IOR1 / IOR2; //take the ratio of the refractive index of air to the object material

    //We know by Snell's Law that sin(theta1)/sin(theta2) = refrRatio
    //We can derive: Rr = (-refrRatio*(normal*ray) - sqrt(1 - refrRatio*refrRatio*(1-(normal*ray)*(normal*ray))))*normal+refrRatio*ray
    double c1 = -ray * normal;
    double c2 = 1 - refrRatio * refrRatio * (1 - c1 * c1);

    if (c2 > 0.0) //if the ray is not at the critical angle to the surface
    {
        c2 = sqrt(c2);
        return refrRatio * ray + (refrRatio * c1 - c2) * normal;
    }

    else
        return NULL;
}

/** 
 * clamp
 *
 * Clamps the given vector to the range [0,1].  Used for colors so that blown out areas or negative color values will not
 * result in an unexpected color.
 *
 * @param  v The vector to be clamped\
 * 
 * @return  The clamped vector
 *
 */
vec3 RayTracer::clamp(vec3 v)
{
    //If any value in the vector is above 1, set it to 1.  If any value is below 0, set it to 0.
    if (v[0] > 1)
        v[0] = 1;
    else if (v[0] < 0)
        v[0] = 0;
    if (v[1] > 1)
        v[1] = 1;
    else if (v[1] < 0)
        v[1] = 0;
    if (v[2] > 1)
        v[2] = 1;
    else if (v[2] < 0)
        v[2] = 0;

    return v;
}

/******************************************************************/
//DESTRUCTOR

RayTracer::~RayTracer()
{
    //Delete all objects in the object list and all lights in the light list
    for (unsigned int i = 0; i < objects.size(); i++)
    {
        delete objects[i];
    }

    for (unsigned int i = 0; i < lights.size(); i++)
    {
        delete lights[i];
    }
}