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

	LightNode cluster(LightNode light1, LightNode light2)
	{
		
	}

	LightNode(PointLight leafLight)
	{
		mainLight = leafLight;
		minBoundBox = new Point(); //Set to 0,0,0
		maxBoundBox = new Point(); // Set to 0,0,0
		isLeaf = true;
	}
};