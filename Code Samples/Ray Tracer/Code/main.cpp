/****************************************************************
 * @file                                                            
 * @author  Samantha Raja <ksam@seas.upenn.edu>                                 
 * @date  Summer 2010
 *
 * @section DESCRIPTION
 * This is the main entry point that creates a RayTracer object, sets the scene, 
 * and triggers the start of the ray tracing.
 *
 ****************************************************************/

#include "raytracer.h"

int main(int argc, char **argv){
    RayTracer r(500,500, 3);
    r.setDefaultScene();
    r.writeImage();
}