#ifndef PBRT_CORE_INTEGRATOR_H
#define PBRT_CORE_INTEGRATOR_H

#include "pbrt.h"
#include "primitive.h"
#include "spectrum.h"
#include "light.h"
#include "reflection.h"
#include "sampler.h"
#include "material.h"
#include "probes.h"
#include "renderer.h"

class LightNode
{
	Light *mainLight;
    int index;
	Spectrum Intensity;
	Point minBoundBox;
	Point maxBoundBox;
	Light *leftChild;
	Light *rightChild;
	bool isLeaf;

	//This will be called while constructing the tree
	LightNode *cluster(LightNode *light1, LightNode *light2);
	
	LightNode(PointLight *leafLight);
};
