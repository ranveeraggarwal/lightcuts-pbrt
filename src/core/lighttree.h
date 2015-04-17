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
	PointLight *mainLight; // Pointer to the representative light
	Spectrum Intensity; // Total intensity in this cluster
	Point minBoundBox; // Bounding Box Co-ordinates
	Point maxBoundBox;
	LightNode *leftChild; // Pointers to children
	LightNode *rightChild;
	bool isLeaf; // Whether this node is a leaf, i.e., reprents a single light

	
	LightNode(PointLight *leafLight); // Leaf constructor
	LightNode(); // Empty constructor
};
LightNode * Cluster(LightNode *light1, LightNode *light2); // Function to create a parent LightNode by clustering two LightNodes
#endif
