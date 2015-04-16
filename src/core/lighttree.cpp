#include "lighttree.h"

LightNode * Cluster(LightNode *light1, LightNode *light2)
{
	LightNode *aCluster = new LightNode();
	//Calculate children
	aCluster->leftChild = light1;
	aCluster->rightChild = light2;
	//New Intensity
	aCluster->Intensity = light1->Intensity + light2->Intensity;
	//The new representative light source 
	aCluster->mainLight = light1->mainLight;
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

LightNode::LightNode(PointLight *leafLight)
	{
		mainLight = leafLight;
		minBoundBox = leafLight->lightPos;
		maxBoundBox = leafLight->lightPos;
        Intensity = leafLight->Intensity;
		isLeaf = true;
	}

LightNode::LightNode() {
}
