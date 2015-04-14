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
	Spectrum Intensity;
	Point minBoundBox;
	Point maxBoundBox;
	Light *leftChild;
	Light *rightChild;
	bool isLeaf;

	//This will be called while constructing the tree
	LightNode cluster(LightNode *light1, LightNode *light2)
	{
		Lightnode *aCluster = new Lightnode();
		//Calculate children
		aCluster->leftChild = light1;
		aCluster->rightChild = light2;
		//New Intensity
		aCluster->Intensity = light1->Intensity + light2.Intensity;
		//The new representative light source 
		if (light1->Intensity > light2->Intensity)
		{
			aCluster->mainLight = light1;
		}
		else
		{
			aCluster->mainLight = light2;
		}
		//It's not a leaf
		aCluster->isLeaf = false;
		//Recalculate the bounding box
		aCluster->minBoundBox.x = min(light1->minBoundBox.x, light2->minBoundBox.x);
		aCluster->minBoundBox.y = min(light1->minBoundBox.y, light2->minBoundBox.y);
		aCluster->minBoundBox.z = min(light1->minBoundBox.z, light2->minBoundBox.z);
		aCluster->maxBoundBox.x = max(light1->maxBoundBox.x, light2->maxBoundBox.x);
		aCluster->maxBoundBox.y = max(light1->maxBoundBox.y, light2->maxBoundBox.y);
		aCluster->maxBoundBox.z = max(light1->maxBoundBox.z, light2->maxBoundBox.z);
		return aCluster;
	}

	LightNode(PointLight *leafLight)
	{
		mainLight = leafLight;
		minBoundBox = leafLight->lightPos;
		maxBoundBox = leafLight->lightPos;
		isLeaf = true;
	}
};