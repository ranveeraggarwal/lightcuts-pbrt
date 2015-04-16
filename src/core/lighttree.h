#ifndef PBRT_CORE_LIGHTTREE_H
#define PBRT_CORE_LIGHTTREE_H

#include "pbrt.h"
#include "primitive.h"
#include "spectrum.h"
#include "light.h"
#include "reflection.h"
#include "sampler.h"
#include "material.h"
#include "probes.h"
#include "renderer.h"
#include "../lights/point.h"

class LightNode
{
public:
	PointLight *mainLight;
    int index;
	Spectrum Intensity;
	Point minBoundBox;
	Point maxBoundBox;
	LightNode *leftChild;
	LightNode *rightChild;
	bool isLeaf;

	//This will be called while constructing the tree
	
	LightNode(PointLight *leafLight);
	LightNode();
};
LightNode * Cluster(LightNode *light1, LightNode *light2);
#endif
