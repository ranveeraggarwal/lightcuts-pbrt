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
	PointLight mainLight;
	Point minBoundBox;
	Point maxBoundBox;
	bool isLeaf;

	//This will be called while constructing the tree
	LightNode cluster(LightNode light1, LightNode light2)
	{
		Lightnode aCluster = new Lightnode();
		if (light1.Intensity > light2.Intensity)
		{
			aCluster.mainLight = light1;
		}
		else
		{
			aCluster.mainLight = light2;
		}
		aCluster.isLeaf = false;
		int dist1 = sqrt((light1.Point.x*light1.Point.x) + (light1.Point.y*light1.Point.y) + (light1.Point.z*light1.Point.z));
		int dist2 = sqrt((light2.Point.x*light2.Point.x) + (light2.Point.y*light2.Point.y) + (light2.Point.z*light2.Point.z));
		if (dist1 < dist2) 
		{
			aCluster.minBoundBox = light1.Point;
			aCluster.maxBoundBox = light2.Point;
		}
		else
		{
			aCluster.minBoundBox = light2.Point;
			aCluster.maxBoundBox = light1.Point;	
		}
		return aCluster;
	}

	LightNode(PointLight leafLight)
	{
		mainLight = leafLight;
		minBoundBox = new Point(); //Set to 0,0,0
		maxBoundBox = new Point(); // Set to 0,0,0
		isLeaf = true;
	}
};