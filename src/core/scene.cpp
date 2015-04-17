
/*
    pbrt source code Copyright(c) 1998-2012 Matt Pharr and Greg Humphreys.

    This file is part of pbrt.

    Redistribution and use in source and binary forms, with or without
    modification, are permitted provided that the following conditions are
    met:

    - Redistributions of source code must retain the above copyright
      notice, this list of conditions and the following disclaimer.

    - Redistributions in binary form must reproduce the above copyright
      notice, this list of conditions and the following disclaimer in the
      documentation and/or other materials provided with the distribution.

    THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS
    IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED
    TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A
    PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
    HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
    SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
    LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
    DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
    THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
    (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
    OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

 */


// core/scene.cpp*
#include "stdafx.h"
#include "scene.h"
#include "camera.h"
#include "film.h"
#include "sampler.h"
#include "volume.h"
#include "parallel.h"
#include "progressreporter.h"
#include "renderer.h"
#include <iostream>
#include <cstdlib>
#include <list>
#include <vector>
using namespace std;
// Scene Method Definitions
Scene::~Scene() {
    delete aggregate;
    delete volumeRegion;
    for (uint32_t i = 0; i < lights.size(); ++i)
        delete lights[i];
}


Scene::Scene(Primitive *accel, const vector<Light *> &lts,
             VolumeRegion *vr) {
    lights = lts;
    aggregate = accel;
    volumeRegion = vr;
    // Scene Constructor Implementation
    bound = aggregate->WorldBound();
    if (volumeRegion) bound = Union(bound, volumeRegion->WorldBound());
}


const BBox &Scene::WorldBound() const {
    return bound;
}

Spectrum CalcCost(LightNode * light1, LightNode * light2) {
    Point minBoundBox;
    Point maxBoundBox;
    // Compute the new bounding box
	minBoundBox.x = min(light1->minBoundBox.x, light2->minBoundBox.x);
	minBoundBox.y = min(light1->minBoundBox.y, light2->minBoundBox.y);
	minBoundBox.z = min(light1->minBoundBox.z, light2->minBoundBox.z);
	maxBoundBox.x = max(light1->maxBoundBox.x, light2->maxBoundBox.x);
	maxBoundBox.y = max(light1->maxBoundBox.y, light2->maxBoundBox.y);
	maxBoundBox.z = max(light1->maxBoundBox.z, light2->maxBoundBox.z);
    double dist;
    dist = (minBoundBox.x - maxBoundBox.x) * (minBoundBox.x - maxBoundBox.x);
    dist += (minBoundBox.y - maxBoundBox.y) * (minBoundBox.y - maxBoundBox.y);
    dist += (minBoundBox.z - maxBoundBox.z) * (minBoundBox.z - maxBoundBox.z);
    return dist * (light1->Intensity + light2->Intensity); // Cost is Sum of intensities multiplied by the square of diagonal length of the bounding box
}

float CompareSpectra(Spectrum &a, Spectrum &b) {
    float x[3];
    a.ToRGB(x);
    float y[3];
    b.ToRGB(y);
    return (x[0] + x[1] + x[2]) - (y[0] + y[1] + y[2]);
}

void Scene::MakeLightTree() {
    cout << "Constructing Light Tree...\n";
    list<int> valid; // List of currently representative lights (root of each seperate subtree). Reduces by 1 on every clustering.
    int n = lights.size();
    vector<LightNode *> v(n); // Vector of upper most LightNode represented by a light
    vector< vector<Spectrum> > costs(n); // Cost matrix for clustering of two LightNodes
    int i;
    for(i = 0; i < n; i++) {
        v[i] = new LightNode((PointLight *)lights[i]); // Initialise with leaf LightNodes
        valid.push_back(i); // Initially, all lights are representative
    }
    int j;
    for(i = 0; i < n; i++) {
        for(j = 0; j < n; j++) {
            Spectrum s(0.);
            if(i >= j) { // To avoid duplication, the matrix is kept upper triangular
                costs[i].push_back(s);
                continue;
            }
            costs[i].push_back(CalcCost(v[i], v[j])); // Initialise cost matrix with pairwise costs
        }
    }
    list<int>::iterator it1;
    while(valid.size() > 1) { // Do until only the root node is left
        int min_i = -1;
        int min_j = -1;
        list<int>::iterator min_iteratori;
        list<int>::iterator min_iteratorj;
        Spectrum min_cost;
        // Find the pair of subtrees with the minimum cost
        for(it1 = valid.begin(); it1 != valid.end(); it1++) {
            list<int>::iterator it2 = it1;
            for(it2++; it2 != valid.end(); it2++) {
                if(min_i == -1 || (CompareSpectra(costs[*it1][*it2], costs[min_i][min_j]) < 0)) {
                    min_i = *it1;
                    min_j = *it2;
                    min_iteratori = it1;
                    min_iteratorj = it2;
                }
            }
        }
        float rgbi[3];
        v[min_i]->Intensity.ToRGB(rgbi);
        float rgbj[3];
        v[min_j]->Intensity.ToRGB(rgbj);
        float pi = rgbi[0] + rgbi[1] + rgbi[2];
        float pj = rgbj[0] + rgbj[1] + rgbj[2];
        int x = 1000 * (pi / (pi + pj));
        // Randomly choose the representative light with probability proportional to the intensity
        if((rand() % 1000) > x) {
            int _temp_ = min_i;
            min_i = min_j;
            min_j = _temp_;
            min_iteratorj = min_iteratori;
        }
        v[min_i] = Cluster(v[min_i], v[min_j]); // Cluster the LightNodes to create a new LightNode
        valid.erase(min_iteratorj); // Remove the non-representative light from the valid list
        for(it1 = valid.begin(); it1 != valid.end(); it1++) { // Update the cost of all affected pairs (The ones with the representative light in the pair)
            if(min_i < *it1) {
                costs[min_i][*it1] = CalcCost(v[min_i], v[*it1]);
            }
            else if(min_i > *it1) {
                costs[*it1][min_i] = CalcCost(v[min_i], v[*it1]);
            }
        }
    }
    it1 = valid.begin();
    lighttree = v[*it1]; // Set the lighttree variable of the Scene object to the root node
    cout << "Light Tree Constructed.\n";
}
